
function rngInit(generator, N)
    #   skipping the expected number is suggested for Sobol
    #   but has been causing spikes for higher particle counts
    if (generator == "Sobol")
        rng = SobolSeq(2)
    elseif (generator == "Random")
        rng = rand(N,2)
    elseif (generator == "Golden")
        rng = GoldenSequence(2)
    else
        print("RNG must be 'Sobol' or 'Random'. ")
    end
    return rng
end

function nextRN(rng, i, generator)
    if (generator == "Sobol")
        tmp_rnd = next!(rng)
        randX = tmp_rnd[1]
        randMu = tmp_rnd[2]
    elseif (generator == "Random")
        randX = rng[i,1]
        randMu = rng[i,2]
    elseif (generator == "Golden")
        randX = rng[i][1]
        randMu = rng[i][2]
    end
    return (randX, randMu)
end


function move_part(midpoints,dx,mu,zone,x,Nx,high_edges,low_edges,dxs,weight,ds_zone,
                    phi_avg, dphi, phi_edge, phi_s, J_avg, J_edge,sigt,
                    exit_right_bins,exit_left_bins,c)
    # direction to sweep
    if mu > 0
        zoneRange = (zone):1:Nx
    else
        zoneRange = (zone):-1:1
    end
    # sweep through all zones
    for z_prop in zoneRange
        # first zone will only be a partial zone length
        if (z_prop == zone)
            if (mu >= 0)
                ds = (high_edges[zone] - x)/mu
            else
                ds = (low_edges[zone] - x)/mu
            end
        # remaining zones will be full length
        else
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

        if (mu > 0)
            J_edge[z_prop+1,:]   += weight
            phi_edge[z_prop+1,:] += weight./abs(mu)
        else
            J_edge[z_prop,:]   -= weight
            phi_edge[z_prop,:] += weight./abs(mu)
        end
    end
        #add to the exiting flux - find the bin and add it to that one
        #exit_left_bins[argmin(abs.(exit_left_bins[:,1] .- mu)),2] += weight/abs(mu)
    return
end
