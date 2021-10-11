include("misc_functions.jl")

function move_part(midpoints,dx,mu,zone,ds_zone,x,Nx,high_edges,low_edges,dxs,weight,
                    phi_avg, dphi, phi_edge, phi_s, J_avg, J_edge,sigt,
                    exit_right_bins,exit_left_bins,c,phi,muSin,z,y,Geo)
    # direction to sweep
    if (muSin*cos(phi)+muSin*sin(phi)+mu > 0)
        zoneRange = (zone):1:Nx
    else
        zoneRange = (zone):-1:1
    end
    # sweep through all zones
    for z_prop in zoneRange
        if (z_prop == zone)
            # first zone will only be a partial zone length
            ds = distance_to_edge(Geo,x,y,z,mu,muSin,phi,low_edges,high_edges,zone)
        else
            # remaining zones will be full length
            ds = ds_zone
        end
        # update variables
        score_TL           = weight.*(1 .- exp.(-(ds*sigt[z_prop,:])))./sigt[z_prop,:]./dxs[z_prop] #implicit capture for track-length
        phi_avg[z_prop,:] += score_TL
        J_avg[z_prop,:]   += score_TL*mu
        phi_s[z_prop,:]   += (weight.*exp.(-low_edges[z_prop]./c)).*(1 .- exp.(-(ds*(sigt[z_prop,:] .+ mu./c))))./(sigt[z_prop,:] .+ mu./c)./dxs[z_prop]
        dphi[z_prop,:]    += (weight.*((-dx/2).*(1 .- exp.(-sigt[z_prop,:]*ds))
                             .+ mu*((1 .- exp.(-sigt[z_prop,:]*ds))./sigt[z_prop,:])
                             .- ds*exp.(-sigt[z_prop,:]*ds)))./sigt[z_prop,:]./dxs[z_prop]
        weight            .*= exp.(-(ds*sigt[z_prop,:]))

        if (muSin*cos(phi)+muSin*sin(phi)+mu > 0)
            J_edge[z_prop+1,:]   += weight
            phi_edge[z_prop+1,:] += weight./abs(mu)
        else
            J_edge[z_prop,:]   -= weight
            phi_edge[z_prop,:] += weight./abs(mu)
        end
        x,y,z = updatePosition(Geo,x,y,z,mu,muSin,phi,ds)
    end
        #add to the exiting flux - find the bin and add it to that one
        #exit_left_bins[argmin(abs.(exit_left_bins[:,1] .- mu)),2] += weight/abs(mu)
    return
end
