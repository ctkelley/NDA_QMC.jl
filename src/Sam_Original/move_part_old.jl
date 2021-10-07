function move_part(midpoints,dx,mu,zone,x,Nx,high_edges,low_edges,dxs,weight,ds_zone,
                     phi_avg, dphi, phi_edge, phi_s, J_avg, J_edge,sigt,
                     exit_right_bins,exit_left_bins,c)
    G = size(sigt)[1]
    if mu > 0
         #get distance to nearest edge
         ds = (high_edges[zone] - x)/abs(mu) #how far in the first zone
         if (sigt[zone,:] > (1e-12)*ones(G))
             score_TL       = weight.*(1 .- exp.(-(ds*sigt[zone,:])))./sigt[zone,:]./dxs[zone,:] #implicit capture for track-length
             phi_avg[zone,:] += score_TL
             J_avg[zone,:]   += score_TL*mu
             phi_s[zone,:]   += (weight.*exp.(-x./c)).*(1 .- exp.(-(ds*(sigt[zone,:] .+ mu./c))))./(sigt[zone,:] .+ mu./c)./dxs[zone,:] #weighted average version

             dphi[zone,:]    += (weight.*((x-midpoints[zone]).*(1 .- exp.(-sigt[zone,:]*ds))
                                .+ mu*((1 .- exp.(-sigt[zone,:]*ds))./sigt[zone,:])
                                .- ds*exp.(-sigt[zone,:]*ds)))./sigt[zone,:]./dxs[zone]
         else
             score_TL       = weight*ds./dxs[zone] #implicit capture for track-length
             phi_avg[zone,:] += score_TL
             J_avg[zone,:]   += score_TL*mu
             phi_s[zone,:]   += phi_avg[zone,:]
         end
         weight           .*= exp.(-(ds*sigt[zone,:])) #decrease the weight
         J_edge[zone+1,:]   += weight
         phi_edge[zone+1,:] += weight./abs(mu)

         #now do all the other zones
         for z_prop in (zone+1):Nx
            if (sigt[z_prop,:] > (1e-12)*ones(G))
                ds += ds_zone
                score_TL           = weight.*(1 .- exp.(-(ds_zone*sigt[z_prop,:])))./sigt[z_prop,:]./dxs[z_prop] #implicit capture for track-length
                phi_avg[z_prop,:] += score_TL
                J_avg[z_prop,:]   += score_TL*mu
                phi_s[z_prop,:]   += (weight.*exp.(-low_edges[z_prop]./c)).*(1 .- exp.(-(ds_zone*(sigt[z_prop,:] .+ mu./c))))./(sigt[z_prop,:] .+ mu./c)./dxs[z_prop]

                dphi[z_prop,:]    += (weight.*((-dx/2).*(1 .- exp.(-sigt[z_prop,:]*ds_zone))
                                     .+ mu*((1 .- exp.(-sigt[z_prop,:]*ds_zone))./sigt[z_prop,:])
                                     .- ds_zone*exp.(-sigt[z_prop,:]*ds_zone)))./sigt[z_prop,:]./dxs[z_prop]
             else
                 score_TL         = weight*ds_zone./dxs[z_prop] #implicit capture for track-length
                 phi_avg[z_prop,:] += score_TL
                 J_avg[z_prop,:]   += score_TL*mu
                 phi_s[z_prop,:]   += phi_avg[z_prop,:]
             end
             weight             .*= exp.(-(ds_zone*sigt[z_prop,:]))
             J_edge[z_prop+1,:]   += weight
             phi_edge[z_prop+1,:] += weight./abs(mu)
         end
         #add to the exiting flux - find the bin and add it to that one
         #exit_right_bins[argmin(abs.(exit_right_bins[:,1] .- mu)),2] += weight/mu
     else
         #get distance to nearest edge
         ds = (low_edges[zone] - x)/mu
         #phi[zone] += weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone]
         if (sigt[zone,:] > (1e-12)*ones(G))
             score_TL       = weight.*(1 .- exp.(-(ds*sigt[zone,:])))./sigt[zone,:]./dxs[zone,:] #implicit capture for track-length
             phi_avg[zone,:] += score_TL
             J_avg[zone,:]   += score_TL*mu
             phi_s[zone,:]   += (weight.*exp.(-x./c)).*(1 .- exp.(-(ds*(sigt[zone,:] .+ mu./c))))./(sigt[zone,:] .+ mu./c)./dxs[zone,:]

             dphi[zone,:]    += (weight.*((x-midpoints[zone]).*(1 .- exp.(-sigt[zone,:]*ds))
                                .+ mu*((1 .- exp.(-sigt[zone,:]*ds))./sigt[zone,:])
                                .- ds*exp.(-sigt[zone,:]*ds)))./sigt[zone,:]./dxs[zone]
         else
             score_TL       = weight*ds./dxs[zone] #implicit capture for track-length
             phi_avg[zone,:] += score_TL
             J_avg[zone,:]   += score_TL*mu
             phi_s[zone,:]   += phi_avg[zone,:]
         end
         weight         .*= exp.(-(ds*sigt[zone,:]))
         J_edge[zone,:]   -= weight
         phi_edge[zone,:] += weight./abs(mu)

         #now do all the other zones
         for z_prop in (zone-1):-1:1
            #phi[z_prop] += weight*(1 - exp(-(ds_zone*sigt[zone])))/sigt[zone]/dxs[zone]
             if (sigt[z_prop,:] > (1e-12)*ones(G))
                 score_TL         = weight.*(1 .- exp.(-(ds_zone*sigt[z_prop,:])))./sigt[z_prop,:]./dxs[z_prop] #implicit capture for track-length
                 phi_avg[z_prop,:] += score_TL
                 J_avg[z_prop,:]   += score_TL*mu
                 phi_s[z_prop,:]   += (weight.*exp.(-high_edges[z_prop]./c)).*(1 .- exp.(-(ds_zone*(sigt[z_prop,:] .+ mu./c))))./(sigt[z_prop,:] .+ mu./c)./dxs[z_prop]

                 dphi[z_prop,:]    += (weight.*((dx/2).*(1 .- exp.(-sigt[z_prop,:]*ds_zone))
                                      .+ mu*((1 .- exp.(-sigt[z_prop,:]*ds_zone))./sigt[z_prop,:])
                                      .- ds_zone*exp.(-sigt[z_prop,:]*ds_zone)))./sigt[z_prop,:]./dxs[z_prop]
             else
                 score_TL         = weight*ds_zone./dxs[z_prop] #implicit capture for track-length
                 phi_avg[z_prop,:] += score_TL
                 J_avg[z_prop,:]   += score_TL*mu
                 phi_s[z_prop,:]   += phi_avg[z_prop,:]
             end
             weight           .*= exp.(-(ds_zone*(sigt[z_prop,:])))
             J_edge[z_prop,:]   -= weight
             phi_edge[z_prop,:] += weight./abs(mu)
         end
         #add to the exiting flux - find the bin and add it to that one
         #exit_left_bins[argmin(abs.(exit_left_bins[:,1] .- mu)),2] += weight/abs(mu)
     end
     return
 end
