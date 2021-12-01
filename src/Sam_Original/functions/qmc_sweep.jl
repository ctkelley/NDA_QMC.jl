
"""
    qmc_sweep(phi_avg, qmc_data)
Performs one sweep of N particles per source (left boundary, right boundary,
volumetric). Returns flux values for iteration.

...
# Arguments
* `phi_avg::Array{Float64,2}`: (Nx)xG matrix of scalar flux values.
* `qmc_data::NamedTuple`: various problem initializiations in qmc_init.jl.
...
"""
function qmc_sweep(phi_avg, qmc_data)

    N = qmc_data.N
    Nx = qmc_data.Nx
    RB = qmc_data.RB
    LB = qmc_data.LB
    G = qmc_data.G
    Geo = qmc_data.Geo

    hasLeft = qmc_data.hasLeft
    hasRight = qmc_data.hasRight

    low_edges = qmc_data.low_edges
    high_edges = qmc_data.high_edges
    midpoints = qmc_data.midpoints
    edges = qmc_data.edges

    dmu = qmc_data.dmu

    exit_left_bins = qmc_data.exit_left_bins
    exit_right_bins = qmc_data.exit_right_bins
    exit_left_bins[:,2] .= 0
    exit_right_bins[:,2] .= 0

    #phi_edge = qmc_data.phi_edge
    #dphi = qmc_data.dphi
    #phi_s = qmc_data.phi_s
    #J_avg = qmc_data.J_avg
    #J_edge = qmc_data.J_edge

    sigt = qmc_data.sigt
    #siga = qmc_data.siga
    source = qmc_data.source
    sigs = qmc_data.sigs
    c = qmc_data.c
    generator =qmc_data.generator

    phi_right = qmc_data.phi_right
    phi_left = qmc_data.phi_left

    #q = phi_avg*sigs' + source # multi group data
    q =  phi_avg.*sigs + source # garcia tests and infinite medium problems
    phi_avg = zeros(Nx,G)
    phi_edge = zeros(Nx+1,G)
    dphi = zeros(Nx,G)
    phi_s = zeros(Nx,G)
    J_avg = zeros(Nx,G)
    J_edge = zeros(Nx+1,G)

    # initialize random number generator, return NxDim matrix
    totalDim = getDim(Geo, hasLeft, hasRight)
    rng = rngInit(generator, N, totalDim)
    # Dim is used to index the rng matrix
    Dim = 1

    if (hasLeft)
        #do left boundary source
        for i in 1:N
            randMu, randPhi = nextBoundaryRN(rng, i, generator, Geo, Dim)
            x = LB + 1e-9          # start at left boundary
            fl(mu) =  sqrt(mu)     # isotropic boundary source (point source is uniform)
            mu = fl(randMu)        # mu in 0 to 1
            #mu = inf_Med_BC_Mu(randMu, x, siga)
            phi = randPhi*2*pi
            y = z = 0
            #determine zone (left boundary)
            zone = 1
            #compute weight
            #weight = q[zone,:]/N*cellVolume(Geo, zone, low_edges, high_edges)*Nx
            weight = phi_left./N*surfaceArea(Geo,x).*ones(G)
            #set phi_edge
            phi_edge[zone,:] = phi_left.*ones(G)#*surfaceArea(Geo,x)
            #total path length across zone
            J_edge[zone,:] += weight
            move_part(  midpoints,mu,x,Nx,high_edges,low_edges,weight,
                        phi_avg, dphi, phi_edge, phi_s, J_avg, J_edge,sigt,
                        exit_right_bins,exit_left_bins,c,phi,z,y,Geo)
        end
        # add to Dim so the next columns in the rng matrix are indexed
        if (Geo == 1)
            Dim += 1
        else
            Dim += 2
        end
    end

    if (hasRight)
        #do right boundary source
        for i in 1:N
            randMu, randPhi = nextBoundaryRN(rng, i, generator, Geo, Dim)
            x = RB - 1e-9         # start at right boundary
            fl(mu) =  -sqrt(mu)
            mu = fl(randMu)       # mu in -1 to 0
            #mu = -inf_Med_BC_Mu(randMu, x, siga)
            phi = randPhi*2*pi
            y = z = 0
            #determine zone (right boundary)
            zone = Nx
            #compute weight
            weight = phi_right./N*surfaceArea(Geo,x).*ones(G)
            #set phi_edge
            phi_edge[zone+1,:] = phi_right.*ones(G)#*surfaceArea(Geo,x)
            #total path length across zone
            J_edge[zone+1,:] += weight
            move_part(  midpoints,mu,x,Nx,high_edges,low_edges,weight,
                        phi_avg, dphi, phi_edge, phi_s, J_avg, J_edge,sigt,
                        exit_right_bins,exit_left_bins,c,phi,z,y,Geo)
        end
        if (Geo == 1)
            Dim += 1
        else
            Dim += 2
        end
    end

    #do volumetric source
    for i in 1:N
        #pick starting point and mu
        randX, randMu, randPhi = nextRN(rng, i, generator, Geo, Dim)
        x =  (RB-LB)*randX    # number between 0 and Lx for starting
        mu = 2*randMu-1       # mu in -1 to 1
        phi = randPhi*2*pi
        y = z = 0
        #compute initial weight
        zone = getZone(x,y,z,low_edges,high_edges)
        #if (any(q[zone,:] .> 1e-12))
        weight = q[zone,:]/N*cellVolume(Geo, zone, low_edges, high_edges)*Nx
        #how far does a particle travel when it crosses a zone
        move_part(  midpoints,mu,x,Nx,high_edges,low_edges,weight,
                    phi_avg, dphi, phi_edge, phi_s, J_avg, J_edge,sigt,
                    exit_right_bins,exit_left_bins,c,phi,z,y,Geo)
        #end
    end

    exit_left_bins[:,2] /= dmu
    exit_right_bins[:,2] /= dmu

    return (phi_avg = phi_avg,
            phi_edge = phi_edge,
            dphi = dphi,
            J_avg = J_avg,
            J_edge = J_edge,
            exit_right_bins = exit_right_bins,
            exit_left_bins = exit_left_bins)
end
