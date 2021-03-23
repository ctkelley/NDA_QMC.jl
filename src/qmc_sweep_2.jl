

function qmc_sweep!(phi_avg, qmc_data)


    function move_part(mu,zone,x,Nx,dt,t,v,high_edges,low_edges,dxs,weight,ds_zone,
    phi_avg, phi_edge, phi_s, J_avg, J_edge,sigt,exit_right_bins, exit_left_bins, c)
       #moving right
       if mu > 0
            # Distance to nearest edge
            ds = (high_edges[zone] - x)/mu #how far in the first zone
            if dt < Inf
                s_max = v*(dt-t) #maximum travel distance given velocity and time increment
            else
                s_max = (Lx - x)/mu
            end
            s = 0
            if (ds >= s_max)
                ds = s_max
            end
            # First partial zone
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
            s += ds

            # Full zones
            z_prop = zone+1
            while (s+ds_zone <= s_max) && (z_prop <= Nx)
                if (sigt[z_prop] > 1e-12)
                    ds += ds_zone
                    score_TL         = weight*(1 - exp(-(ds_zone*sigt[z_prop])))/sigt[z_prop]/dxs[z_prop] #implicit capture for track-length
                    phi_avg[z_prop] += score_TL
                    J_avg[z_prop]   += score_TL*mu
                    phi_s[z_prop] += (weight*exp(-low_edges[z_prop]/c))*(1 - exp(-(ds_zone*(sigt[z_prop] + mu/c))))/(sigt[z_prop] + mu/c)/dxs[z_prop]
                    s += ds_zone
                else
                    score_TL         = weight*ds_zone/dxs[z_prop] #implicit capture for track-length
                    phi_avg[z_prop] += score_TL
                    J_avg[z_prop]   += score_TL*mu
                    phi_s[z_prop] += phi_avg[z_prop]
                    s += ds_zone
                end
                weight             *= exp(-(ds_zone*sigt[z_prop]))
                J_edge[z_prop+1]   += weight
                phi_edge[z_prop+1] += weight/abs(mu)
                z_prop += 1
            end

            # Last partial zone
            if (z_prop < Nx) && (s < s_max)
                zone = argmax(1*((s*mu).>=low_edges).*((s*mu) .< high_edges))
                ds = s_max-s
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
                J_edge[zone+1]   += weight #current density tally
                phi_edge[zone+1] += weight/abs(mu)
            end
            #add to the exiting flux - find the bin and add it to that one
            exit_right_bins[argmin(abs.(exit_right_bins[:,1] .- mu)),2] += weight/mu
        #moving left
        else
            #get distance to nearest edge
            ds = (low_edges[zone] - x)/mu
            if dt < Inf
                s_max = v*(dt-t) #maximum travel distance given velocity and time increment
            else
                s_max = -x/mu + 1e-6
            end
            s = 0
            if (ds >= s_max)
                ds = s_max
            end
            # First partial Zone
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
            s += ds
            # Full zones
            z_prop = zone-1
            while (s+ds_zone <= s_max) && (z_prop >= 1)
                if (sigt[z_prop] > 1e-12)
                    score_TL         = weight*(1 - exp(-(ds_zone*sigt[z_prop])))/sigt[z_prop]/dxs[z_prop] #implicit capture for track-length
                    phi_avg[z_prop] += score_TL
                    J_avg[z_prop]   += score_TL*mu
                    phi_s[z_prop] += (weight*exp(-high_edges[z_prop]/c))*(1 - exp(-(ds_zone*(sigt[z_prop] + abs(mu)/c))))/(sigt[z_prop] + abs(mu)/c)/dxs[z_prop]
                    s += ds_zone
                else
                    score_TL         = weight*ds_zone/dxs[z_prop] #implicit capture for track-length
                    phi_avg[z_prop] += score_TL
                    J_avg[z_prop]   += score_TL*mu
                    phi_s[z_prop] += phi_avg[z_prop]
                    s += ds_zone
                end
               weight           *= exp(-(ds_zone*(sigt[z_prop])))
               J_edge[z_prop]   -= weight
               phi_edge[z_prop] += weight/abs(mu)
               z_prop -= 1
           end

            # Last partial zone
            if (z_prop > 1) && (s < s_max)
                zone = argmax(1*((s*mu).>=low_edges).*((s*mu) .< high_edges))
                ds = (low_edges[zone] - s*mu)/mu
                if (sigt[zone] > 1e-12)
                    score_TL       = weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone] #implicit capture for track-length
                    phi_avg[zone] += score_TL
                    J_avg[zone]   += score_TL*mu
                    phi_s[zone] += (weight*exp(-x/c))*(1 - exp(-(ds*(sigt[zone] + abs(mu)/c))))/(sigt[zone] + abs(mu)/c)/dxs[zone]
                else
                    score_TL       = weight*ds/dxs[zone] #implicit capture for track-length
                    phi_avg[zone] += score_TL
                    J_avg[zone]   += score_TL*mu
                    phi_s[zone] += phi_edge[zone]
                end
                weight         *= exp.(-(ds*sigt[zone]))
                J_edge[zone]   -= weight
                phi_edge[zone] += weight/abs(mu)
            end
            #add to the exiting flux - find the bin and add it to that one
            exit_left_bins[argmin(abs.(exit_left_bins[:,1] .- mu)),2] += weight/abs(mu)
        end
        if dt < Inf
            print(true)
            phi_avg /= dt
            phi_edge /= dt
            J_avg /= dt
            J_edge /= dt
            phi_s /= dt
        end
        return (weight, phi_avg, phi_edge, phi_s, J_avg, J_edge, exit_right_bins,
                                                                    exit_left_bins)
    end

    #do boundary source
    for i in 1:N
        tmp_rnd = next!(rng_bndl) #for Sobol
        x = 0 #start at left boundary
        mu = fl(tmp_rnd[1])  #mu in 0 to 1
        if dt == Inf
            t = 0
        else
            t = dt*fl(tmp_rnd[2]) # start time somewhere in dt
        end
        #determine zone (left boundary)
        zone = 1
        #compute weight
        weight = 1/Nb*0.5
        #set phi_edge
        phi_edge[1] = 1
        #total path length across zone
        ds_zone = dx/abs(mu)
        J_edge[zone] += weight
        move_part(mu,zone,x,Nx,dt,t,v,high_edges,low_edges,dxs,weight,ds_zone,
                                phi_avg, phi_edge, phi_s, J_avg, J_edge,sigt,
                                        exit_right_bins,exit_left_bins,c)
    end
    
    #do sources in the problem
    for i in 1:N
        #pick starting point and mu
        tmp_rnd = next!(rng) #for Sobol
        x = Lx*tmp_rnd[1] # number between 0 and Lx for starting
        mu = 2*tmp_rnd[2]-1 #mu in -1 to 1
        if dt == Inf
            t = 0
        else
            t = dt*tmp_rnd[3] # start time somewhere in dt
        end
        #determine zone
        zone = argmax(1*(x.>=low_edges).*(x .< high_edges))
        #compute weight
        weight = q[zone]/N*dx*Nx
        #how far does a particle travel when it crosses a zone
        ds_zone = dx/abs(mu)
        move_part(mu,zone,x,Nx,dt,t,v,high_edges,low_edges,dxs,weight,ds_zone,
                                phi_avg, phi_edge, phi_s, J_avg, J_edge,sigt,
                                        exit_right_bins,exit_left_bins,c)
    end

    #first order forward difference
    for i in 1:(Nx-1)
        #spatial derivative of scalar flux
        dphi[i] = phi_edge[i+1] - phi_edge[i]
    end

    dphi /= (dx)
    exit_left_bins[:,2] /= dmu
    exit_right_bins[:,2] /= dmu

    return phi_avg, phi_edge, J_avg, J_edge, psi_right, psi_left
