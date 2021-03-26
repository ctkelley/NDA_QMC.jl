
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
