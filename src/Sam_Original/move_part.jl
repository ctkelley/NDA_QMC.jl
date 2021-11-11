include("misc_functions.jl")

function move_part(midpoints,mu,x,Nx,high_edges,low_edges,weight,
                    phi_avg, dphi, phi_edge, phi_s, J_avg, J_edge,sigt,
                    exit_right_bins,exit_left_bins,c,phi,z,y,Geo)

    # collect list of zones crossed and path lengths in each
    zoneList, dsList = getPath(Geo,x,y,z,mu,phi,low_edges,high_edges)
    counter = 1

    # sweep through all zones
    for z_prop in zoneList
        ds = dsList[counter]
        dV = cellVolume(Geo,z_prop,low_edges,high_edges)
        # update variables
        score_TL           = weight.*(1 .- exp.(-(ds*sigt[z_prop,:])))./(sigt[z_prop,:].*dV) #implicit capture for track-length
        phi_avg[z_prop,:] += score_TL
        J_avg[z_prop,:]   += score_TL*mu
        phi_s[z_prop,:]   += (weight.*exp.(-low_edges[z_prop]./c)).*(1 .- exp.(-(ds*(sigt[z_prop,:] .+ mu./c))))./(sigt[z_prop,:] .+ mu./c)./cellVolume(Geo,z_prop,low_edges,high_edges)
        dphi[z_prop,:]    += (weight.*((-dV/2).*(1 .- exp.(-sigt[z_prop,:]*ds))
                             .+ mu*((1 .- exp.(-sigt[z_prop,:]*ds))./sigt[z_prop,:])
                             .- ds*exp.(-sigt[z_prop,:]*ds)))./sigt[z_prop,:]./cellVolume(Geo,z_prop,low_edges,high_edges)
        weight            .*= exp.(-(ds*sigt[z_prop,:]))

        if (mu > 0)
            J_edge[z_prop+1,:]   += weight
            phi_edge[z_prop+1,:] += weight./abs(mu)
        else
            J_edge[z_prop,:]   -= weight
            phi_edge[z_prop,:] += weight./abs(mu)
        end
        counter += 1
    end
        #add to the exiting flux - find the bin and add it to that one
        #exit_left_bins[argmin(abs.(exit_left_bins[:,1] .- mu)),2] += weight/abs(mu)
    return
end
