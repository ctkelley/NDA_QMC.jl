
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

function distance_to_edge(Geo,x,y,z,mu,muSin,phi,low_edges,high_edges,zone)
    # which edge is the particle traveling to
    ds = 0
    if (mu >= 0)
        edge = high_edges[zone]
    else
        edge = low_edges[zone]
    end
    # distance to edge calc varies by geometry
    if (Geo == 1) # slab
        ds = (edge - x)/mu
    else
        if (Geo == 2) # cylinder
            a = muSin*cos(phi)^2 + muSin*sin(phi)^2
            h = (x*muSin*cos(phi) + y*muSin*sin(phi))
        else # sphere
            a = 1
            h = (z*muSin*cos(phi) + y*muSin*sin(phi) + x*mu)
        end
        c = x^2 + y^2 + z^2 - edge^2
        if ((a != 0) && (h^2 - a*c > 0))
            d1 = (-h + sqrt(h^2 - a*c))/a
            d2 = (-h - sqrt(h^2 - a*c))/a
            if (c < 0)
                if (d1 > 0)
                    ds = d1
                else
                    ds = d2
                end
            elseif (d1 > 0)
                if (d1 > d2)
                    ds = d2
                else
                    ds = d1
                end
            else
                println(" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ")
                println("%%%%%% No intersection %%%%%% ")
                println(x, ", ", y, ", ",z)
                println(mu, ", ", muSin, ", ", phi)
                println(zone)
                println(" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ")

            end
        end
    end
    return ds
end


function distance_across_zone(Geo,x,y,z,mu,muSin,phi,low_edges,high_edges,zone,dxs)
    if (muSin*cos(phi)+muSin*sin(phi)+mu >= 0)
        r2 = high_edges[zone]
        r1 = low_edges[zone]
    else
        r2 = low_edges[zone]
        r1 = high_edges[zone]
    end
    ds = sqrt( (r2-r1)^2 / ((muSin*cos(phi))^2 + (muSin*sin(phi))^2 + mu^2) )
    return ds
end


function updatePosition(Geo,x,y,z,mu,muSin,phi,ds)
    if (Geo == 2)
        #CYLINDER
        x += ds*muSin*cos(phi)
        y += ds*muSin*sin(phi)
    else
        #SLAB
        x += ds*mu
        if (Geo == 3)
            #SPHERE
            z += ds*muSin*cos(phi)
            y += ds*muSin*sin(phi)
        end
    end
    return (x,y,z)
end

function getZone(x,y,z,low_edges,high_edges)
    r = sqrt(x^2 + y^2 + z^2)
    return argmax(1*(r.>=low_edges).*(r .< high_edges))
end
