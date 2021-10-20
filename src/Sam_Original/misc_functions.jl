
function rngInit(generator, Geo, N)
    #   skipping the expected number is suggested for Sobol
    #   but has been causing spikes for higher particle counts
    if (Geo == 1)
        Dim = 2
    else
        Dim = 3
    end

    if (generator == "Sobol")
        rng = SobolSeq(Dim)
    elseif (generator == "Random")
        rng = rand(N,Dim)
    elseif (generator == "Golden")
        rng = GoldenSequence(Dim)
    else
        print("RNG must be 'Sobol', 'Golden', or 'Random'. ")
    end
    return rng
end

function nextRN(rng, i, generator, Geo)
    if (Geo == 1)
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
        randPhi = 0
    else
        if (generator == "Sobol")
            tmp_rnd = next!(rng)
            randX = tmp_rnd[1]
            randMu = tmp_rnd[2]
            randPhi = tmp_rnd[3]
        elseif (generator == "Random")
            randX = rng[i,1]
            randMu = rng[i,2]
            randPhi = rng[i,3]
        elseif (generator == "Golden")
            randX = rng[i][1]
            randMu = rng[i][2]
            randPhi = rng[i][3]
        end
    end
    return (randX, randMu, randPhi)
end

function cellVolume(Geo, zone, low_edges, high_edges)
    if (Geo == 1)
        return (high_edges[zone] - low_edges[zone])
    elseif (Geo == 2)
        return pi*(high_edges[zone]^2 - low_edges[zone]^2)
    elseif (Geo == 3)
        return (4/3)*pi*(high_edges[zone]^3 - low_edges[zone]^3)
    end
end


function sphere_edge(x,y,z,mu,muSin,phi,edge)
    ds = Inf
    #h = (z*muSin*cos(phi) + y*muSin*sin(phi) + x*mu)
    h = (x*cos(phi) + y*sin(phi) + z*mu)
    c = x^2 + y^2 + z^2 - edge^2

    if (h^2 - c > 0)
        d1 = (-h + sqrt(h^2 - c))
        d2 = (-h - sqrt(h^2 - c))
        if (c < 0)
            # "inside" the sphere
            if (d1 > 0)
                ds = d1
            else
                ds = d2
            end
        elseif (d1 > 0)
            if ((d1 > d2) && (d2 > 0))
                ds = d2
            else
                ds = d1
            end
        end
    end
    return ds
end

function cylinder_edge(x,y,mu,muSin,phi,edge) #get rid of mu
    ds = Inf
    #a = muSin*cos(phi)^2 + muSin*sin(phi)^2
    #h = (x*muSin*cos(phi) + y*muSin*sin(phi))
    a = cos(phi)^2 + sin(phi)^2
    h = (x*cos(phi) + y*sin(phi))
    c = x^2 + y^2 - edge^2

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
            if ((d1 > d2) && (d2 > 0))# could d2 be negative
                ds = d2
            else
                ds = d1
            end
        end
    end
    return ds
end

function slab_edge(mu,x,edge)
    return (edge - x)/mu
end

function distance_to_edge(Geo,x,y,z,mu,muSin,phi,low_edges,high_edges,zone)
    # which edge is the particle traveling to
    # distance to edge calc varies by geometry
    if (Geo == 1) # slab
        if (mu >= 0)
            ds = slab_edge(mu,x,high_edges[zone])
        else
            ds = slab_edge(mu,x,low_edges[zone])
        end
    elseif (Geo == 2) # cyliner
        ds1 = cylinder_edge(x,y,mu,muSin,phi,high_edges[zone])
        ds2 = cylinder_edge(x,y,mu,muSin,phi,low_edges[zone])
        if (abs(ds1)<abs(ds2))
            ds = ds1
        elseif (abs(ds1)>abs(ds2))
            ds = ds2
        end
    elseif (Geo == 3)
        ds1 = sphere_edge(x,y,z,mu,muSin,phi,high_edges[zone])
        ds2 = sphere_edge(x,y,z,mu,muSin,phi,low_edges[zone])
        if (abs(ds1)<abs(ds2))
            ds = ds1
        elseif (abs(ds1)>abs(ds2))
            ds = ds2
        end
    end
    return (ds)
end

function updatePosition(Geo,x,y,z,mu,muSin,phi,ds)
    if (Geo == 2)
        #CYLINDER
        x += ds*cos(phi)
        y += ds*sin(phi)
    else
        #SLAB
        z += ds*mu
        if (Geo == 3)
            #SPHERE
            x += ds*cos(phi)
            y += ds*sin(phi)
        end
    end
    return (x,y,z)
end

function getZone(x,y,z,low_edges,high_edges)
    r = sqrt(x^2 + y^2 + z^2)
    return argmax(1*(r.>=low_edges).*(r .< high_edges))
end

function getNextZone(Geo,x,y,z,mu,muSin,phi,ds,low_edges,high_edges)
    x,y,z = updatePosition(Geo,x,y,z,mu,muSin,phi,ds)
    r = sqrt(x^2 + y^2 + z^2)
    return argmax(1*(r.>=low_edges).*(r .< high_edges))
end

function getRadius(x,y,z)
    return sqrt(x^2 + y^2 + z^2)
end

function getPath(Geo,x,y,z,mu,muSin,phi,low_edges,high_edges)
    zone = getZone(x,y,z,low_edges,high_edges)
    if (Geo == 1)
        # direction to sweep
        if (mu >= 0)
            zoneList = (zone):1:Nx # going right
        else
            zoneList = (zone):-1:1 # going left
        end
        # partial zone
        ds = distance_to_edge(Geo,x,y,z,mu,muSin,phi,low_edges,high_edges,zone)
        # full zones
        dsList = ones(size(zoneList)[1])*abs(cellVolume(Geo,zone,low_edges,high_edges)/mu)
        dsList[1] = ds
    else
        alive = true
        dsList = []
        zoneList = []
        # here we dont know exactly which zones it will cross, so calculate one
        # at a time and move the particle
        while (alive)
            zone = getZone(x,y,z,low_edges,high_edges)
            ds = distance_to_edge(Geo,x,y,z,mu,muSin,phi,low_edges,high_edges,zone) + 1e-9
            append!(dsList, ds)
            append!(zoneList, zone)
            x,y,z = updatePosition(Geo,x,y,z,mu,muSin,phi,ds)
            if (zone == Nx)
                alive = false
            end
        end
    end
    return (zoneList, dsList)
end
