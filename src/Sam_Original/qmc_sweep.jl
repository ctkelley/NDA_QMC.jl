include("move_part.jl")

function qmc_sweep(phi_avg, qmc_data)

    N = qmc_data.N
    Nx = qmc_data.Nx
    Lx = qmc_data.Lx
    dx = qmc_data.dx
    G = qmc_data.G

    low_edges = qmc_data.low_edges
    high_edges = qmc_data.high_edges
    dxs = qmc_data.dxs
    midpoints = qmc_data.midpoints
    edges = qmc_data.edges

    dmu = qmc_data.dmu

    exit_left_bins = qmc_data.exit_left_bins
    exit_right_bins = qmc_data.exit_right_bins
    exit_left_bins[:,2] .= 0
    exit_right_bins[:,2] .= 0

    phi_edge = qmc_data.phi_edge
    dphi = qmc_data.dphi
    phi_s = qmc_data.phi_s
    J_avg = qmc_data.J_avg
    J_edge = qmc_data.J_edge

    sigt = qmc_data.sigt
    source = qmc_data.source
    sigs = qmc_data.sigs
    c = qmc_data.c

    q = phi_avg*sigs' + source
    phi_avg = zeros(Nx,G)

    #initialize sobol sequence
    #   skipping the expected number is suggested for Sobol
    #   but has been causing spikes for higher particle counts
    rng = SobolSeq(2)
    rng_bndl = SobolSeq(1)
    #skip(rng,N)
    #skip(rng_bndl,N)

    """
    #do left boundary source
    for i in 1:N
        tmp_rnd = next!(rng_bndl) #for Sobol
        x = 0                     #start at left boundary
        fl(mu) =  sqrt(mu)
        mu = fl(tmp_rnd[1])       #mu in 0 to 1
        #determine zone (left boundary)
        zone = 1
        #compute weight
        weight = ones(G)*0.5/(N)#1/N*0.5
        #set phi_edge
        phi_edge[1,:] = ones(G)#1
        #total path length across zone
        ds_zone = dx/abs(mu)
        J_edge[zone,:] += weight
        move_part(midpoints,dx,mu,zone,x,Nx,high_edges,low_edges,dxs,weight,ds_zone,
                  phi_avg, dphi, phi_edge, phi_s, J_avg, J_edge,sigt,
                  exit_right_bins,exit_left_bins,c)
    end
    """
    """
    #do right boundary source
    for i in 1:N
        tmp_rnd = next!(rng_bndl) #for Sobol
        x = Lx                     #start at left boundary
        fl(mu) =  -sqrt(mu)
        mu = fl(tmp_rnd[1])       #mu in -1 to 0
        #determine zone (left boundary)
        zone = Nx
        #compute weight
        weight = ones(G)*0.5/(N)#1/N*0.5
        #set phi_edge
        phi_edge[Nx,:] = ones(G)#1
        #total path length across zone
        ds_zone = dx/abs(mu)
        J_edge[zone,:] += weight
        move_part(mu,zone,x,Nx,high_edges,low_edges,dxs,weight,ds_zone,
                  phi_avg, dphi, phi_edge, phi_s, J_avg, J_edge,sigt,
                  exit_right_bins,exit_left_bins,c)
    end
    """
    #do volumetric source
    for i in 1:N
        #pick starting point and mu
        tmp_rnd = next!(rng)    # for Sobol
        x = Lx*tmp_rnd[1]       # number between 0 and Lx for starting
        mu = 2*tmp_rnd[2]-1     # mu in -1 to 1
        #determine zone
        zone = argmax(1*(x.>=low_edges).*(x .< high_edges))
        #compute initial weight
        weight = q[zone,:]/N*dx*Nx
        #how far does a particle travel when it crosses a zone
        ds_zone = dx/abs(mu)
        move_part(midpoints,dx,mu,zone,x,Nx,high_edges,low_edges,dxs,weight,ds_zone,
                  phi_avg, dphi, phi_edge, phi_s, J_avg, J_edge,sigt,
                  exit_right_bins,exit_left_bins,c)
    end

    return (phi_avg = phi_avg,
            phi_edge = phi_edge,
            dphi = dphi,
            J_avg = J_avg,
            J_edge = J_edge)
end
