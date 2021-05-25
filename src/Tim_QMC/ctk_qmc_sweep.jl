#include("qmc_move_part.jl")

function qmc_sweep(phi_avg, qmc_data)
    
    function move_part(mu,zone,x,Nx,high_edges,low_edges,dxs,weight,ds_zone,
                                    phi_avg, phi_edge, phi_s, J_avg, J_edge,sigt,
                                            exit_right_bins,exit_left_bins,c)
       if mu > 0
            #get distance to nearest edge
            ds = (high_edges[zone] - x)/abs(mu) #how far in the first zone
            if (sigt[zone] > 1e-12)
                score_TL       = weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone] #implicit capture for track-length
                phi_avg[zone] += score_TL
                J_avg[zone]   += score_TL*mu
                phi_s[zone] += (weight*exp(-x/c))*(1 - exp(-(ds*(sigt[zone] + mu/c))))/(sigt[zone] + mu/c)/dxs[zone] #weighted average version
            else
                score_TL       = weight*ds/dxs[zone] #implicit capture for track-length
                phi_avg[zone] += score_TL
                J_avg[zone]   += score_TL*mu
                phi_s[zone] += phi_avg[zone]
            end
            weight           *= exp(-(ds*sigt[zone])) #decrease the weight
            J_edge[zone+1]   += weight
            phi_edge[zone+1] += weight/abs(mu)

            #now do all the other zones
            for z_prop in (zone+1):Nx
               #phi[z_prop] += weight*(1 - exp(-(ds_zone*sigt[zone])))/sigt[zone]/dxs[zone]
               if (sigt[z_prop] > 1e-12)
                   ds += ds_zone
                   score_TL         = weight*(1 - exp(-(ds_zone*sigt[z_prop])))/sigt[z_prop]/dxs[z_prop] #implicit capture for track-length
                   phi_avg[z_prop] += score_TL
                   J_avg[z_prop]   += score_TL*mu
                   phi_s[z_prop] += (weight*exp(-low_edges[z_prop]/c))*(1 - exp(-(ds_zone*(sigt[z_prop] + mu/c))))/(sigt[z_prop] + mu/c)/dxs[z_prop]
                else
                    score_TL         = weight*ds_zone/dxs[z_prop] #implicit capture for track-length
                    phi_avg[z_prop] += score_TL
                    J_avg[z_prop]   += score_TL*mu
                    phi_s[z_prop] += phi_avg[z_prop]
                end
                weight             *= exp(-(ds_zone*sigt[z_prop]))
                J_edge[z_prop+1]   += weight
                phi_edge[z_prop+1] += weight/abs(mu)
            end

            #add to the exiting flux - find the bin and add it to that one
            exit_right_bins[argmin(abs.(exit_right_bins[:,1] .- mu)),2] += weight/mu
        else
            #get distance to nearest edge
            ds = (low_edges[zone] - x)/mu
            #phi[zone] += weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone]
            if (sigt[zone] > 1e-12)
                score_TL       = weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone] #implicit capture for track-length
                phi_avg[zone] += score_TL
                J_avg[zone]   += score_TL*mu
                phi_s[zone] += (weight*exp(-x/c))*(1 - exp(-(ds*(sigt[zone] + abs(mu)/c))))/(sigt[zone] + abs(mu)/c)/dxs[zone]
            else
                score_TL       = weight*ds/dxs[zone]
                phi_avg[zone] += score_TL
                J_avg[zone]   += score_TL*mu
                phi_s[zone] += phi_avg[zone]
            end
            weight         *= exp.(-(ds*sigt[zone]))
            J_edge[zone]   -= weight
            phi_edge[zone] += weight/abs(mu)

            #now do all the other zones
            for z_prop in (zone-1):-1:1
               #phi[z_prop] += weight*(1 - exp(-(ds_zone*sigt[zone])))/sigt[zone]/dxs[zone]
                if (sigt[z_prop] > 1e-12)
                    score_TL         = weight*(1 - exp(-(ds_zone*sigt[z_prop])))/sigt[z_prop]/dxs[z_prop] #implicit capture for track-length
                    phi_avg[z_prop] += score_TL
                    J_avg[z_prop]   += score_TL*mu
                    phi_s[z_prop] += (weight*exp(-high_edges[z_prop]/c))*(1 - exp(-(ds_zone*(sigt[z_prop] + abs(mu)/c))))/(sigt[z_prop] + abs(mu)/c)/dxs[z_prop]
                else
                    score_TL         = weight*ds_zone/dxs[z_prop] #implicit capture for track-length
                    phi_avg[z_prop] += score_TL
                    J_avg[z_prop]   += score_TL*mu
                    phi_s[z_prop] += phi_avg[z_prop]
                end
                weight           *= exp(-(ds_zone*(sigt[z_prop])))
                J_edge[z_prop]   -= weight
                phi_edge[z_prop] += weight/abs(mu)
            end

            #add to the exiting flux - find the bin and add it to that one
            exit_left_bins[argmin(abs.(exit_left_bins[:,1] .- mu)),2] += weight/abs(mu)
        end
        return
    end

    N = qmc_data.N
    Nx = qmc_data.Nx
    Lx = qmc_data.Lx
    dx = qmc_data.dx

    low_edges = qmc_data.low_edges
    high_edges = qmc_data.high_edges
    dxs = qmc_data.dxs
    midpoints = qmc_data.midpoints
    edges = qmc_data.edges

    #rng = qmc_data.rng
    #rng_bndl = qmc_data.rng_bndl
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
    q = phi_avg.*sigs .+ source
    phi_avg = zeros(Nx)

    #initialize sobol sequence
    rng = SobolSeq(2)
    skip(rng,N) #skipping the expected number is suggested for Sobol
    rng_bndl = SobolSeq(1)
    skip(rng_bndl,N) #skipping the expected number is suggested for Sobol

    #do boundary source
    for i in 1:N
        tmp_rnd = next!(rng_bndl) #for Sobol
        x = 0 #start at left boundary
        fl(mu) =  sqrt(mu)
        mu = fl(tmp_rnd[1])  #mu in 0 to 1
        #determine zone (left boundary)
        zone = 1
        #compute weight
        weight = 1/N*0.5
        #set phi_edge
        phi_edge[1] = 1
        #total path length across zone
        ds_zone = dx/abs(mu)
        J_edge[zone] += weight
        move_part(mu,zone,x,Nx,high_edges,low_edges,dxs,weight,ds_zone,
                                        phi_avg, phi_edge, phi_s, J_avg, J_edge,sigt,
                                                exit_right_bins,exit_left_bins,c)
    end

    #do sources in the problem
    for i in 1:N
        #pick starting point and mu
        tmp_rnd = next!(rng) #for Sobol
        x = Lx*tmp_rnd[1] # number between 0 and Lx for starting
        mu = 2*tmp_rnd[2]-1 #mu in -1 to 1
        #determine zone
        zone = argmax(1*(x.>=low_edges).*(x .< high_edges))
        #compute weight
        weight = q[zone]/N*dx*Nx
        #how far does a particle travel when it crosses a zone
        ds_zone = dx/abs(mu)
        move_part(mu,zone,x,Nx,high_edges,low_edges,dxs,weight,ds_zone,
                                        phi_avg, phi_edge, phi_s, J_avg, J_edge,sigt,
                                                exit_right_bins,exit_left_bins,c)
    end

    # spatial derivative of phi
    for i in 1:(Nx-1)
        dphi[i] = phi_edge[i+1] - phi_edge[i]
    end

    dphi ./= (dx)
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
