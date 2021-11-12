
function rngInit(generator, Geo, N, Dim)
    #   skipping the expected number is suggested for Sobol
    #   but has been causing spikes for higher particle counts

    # need to  create a matrix of random numbers by number of dimensions
    # the sobol package makes this difficult
    if (generator == "Sobol")
        sobol = SobolSeq(Dim)
        rng = zeros(N, Dim)
        for i in 1:N
            temp = next!(sobol)
            rng[i,:] = temp
        end
    elseif (generator == "Random")
        rng = rand(N,Dim)
    elseif (generator == "Golden")
        rng = GoldenSequence(Dim)
    else
        println()
        println("RNG must be 'Sobol', 'Golden', or 'Random'. ")
        println()
    end
    return rng
end
function nextBoundaryRN(rng, i, generator, Geo, Dim)
    if (Geo == 1)
        if (generator == "Sobol")
            randMu = rng[i,Dim] .+ 1e-9
        elseif (generator == "Random")
            randMu = rng[i,Dim]
        elseif (generator == "Golden")
            randMu = rng[i][Dim]
        end
        randPhi = 0
    else
        if (generator == "Sobol")
            randMu = rng[i,Dim]
            randPhi = rng[i,Dim+1]
        elseif (generator == "Random")
            randMu = rng[i,Dim]
            randPhi = rng[i,Dim+1]
        elseif (generator == "Golden")
            randMu = rng[i][Dim]
            randPhi = rng[i][Dim+1]
        end
    end
    return (randMu, randPhi)
end

function nextRN(rng, i, generator, Geo, Dim)
    if (Geo == 1)
        if (generator == "Sobol")
            randX = rng[i,Dim] .+ 1e-9
            randMu = rng[i,Dim+1].+ 1e-9
        elseif (generator == "Random")
            randX = rng[i,Dim]
            randMu = rng[i,Dim+1]
        elseif (generator == "Golden")
            randX = rng[i][Dim]
            randMu = rng[i][Dim+1]
        end
        randPhi = 0
    else
        if (generator == "Sobol")
            randX = rng[i,Dim]
            randMu = rng[i,Dim+1]
            randPhi = rng[i,Dim+2]
        elseif (generator == "Random")
            randX = rng[i,Dim]
            randMu = rng[i,Dim+1]
            randPhi = rng[i,Dim+2]
        elseif (generator == "Golden")
            randX = rng[i][Dim]
            randMu = rng[i][Dim+1]
            randPhi = rng[i][Dim+2]
        end
    end
    return (randX, randMu, randPhi)
end

function cellVolume(Geo, zone, low_edges, high_edges)
    if (Geo == 1) # slab
        return (high_edges[zone] - low_edges[zone])
    elseif (Geo == 2) # cylinder
        return pi*(high_edges[zone]^2 - low_edges[zone]^2)
    elseif (Geo == 3) # sphere
        return (4/3)*pi*(high_edges[zone]^3 - low_edges[zone]^3)
    end
end

function surfaceArea(Geo, r)
    if (Geo == 1) # slab
        return (0.5)
    elseif (Geo == 2) # cylinder
        return 2*pi*(r)
    elseif (Geo == 3) # sphere
        return (4)*pi*(r^2)
    end
end


function sphere_edge(x,y,z,mu,phi,edge)
    ds = Inf
    #h = (x*muSin*cos(phi) + y*muSin*sin(phi) + z*mu)
    h = (x*sqrt(1-mu^2)*cos(phi) + y*sqrt(1-mu^2)*sin(phi) + z*mu)
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

function cylinder_edge(x,y,mu,phi,edge) #get rid of mu
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

function distance_to_edge(Geo,x,y,z,mu,phi,low_edges,high_edges,zone)
    # which edge is the particle traveling to
    # distance to edge calc varies by geometry
    if (Geo == 1) # slab
        if (mu >= 0)
            ds = slab_edge(mu,x,high_edges[zone])
        else
            ds = slab_edge(mu,x,low_edges[zone])
        end
    elseif (Geo == 2) # cyliner
        ds1 = cylinder_edge(x,y,mu,phi,high_edges[zone])
        ds2 = cylinder_edge(x,y,mu,phi,low_edges[zone])
        if (abs(ds1)<=abs(ds2))
            ds = ds1
        elseif (abs(ds1)>abs(ds2))
            ds = ds2
        end
    elseif (Geo == 3)
        ds1 = sphere_edge(x,y,z,mu,phi,high_edges[zone])
        ds2 = sphere_edge(x,y,z,mu,phi,low_edges[zone])
        if (abs(ds1)<=abs(ds2))
            ds = ds1
        elseif (abs(ds1)>abs(ds2))
            ds = ds2
        end
    end
    if (!@isdefined ds)
        println()
        println("ds1 = ", ds1)
        println("ds2 = ", ds2)
        println()
        println("x = ", x)
        println("y = ", y)
        println("z = ", z)
        println("mu = ", mu)
        println("phi = ", phi)
        println("zone = ", zone)
        println()
    end
    return (ds)
end

function updatePosition(Geo,x,y,z,mu,phi,ds)
    if (Geo == 2)
        #CYLINDER
        x += ds*cos(phi)
        y += ds*sin(phi)
    else
        #SLAB
        z += ds*mu
        if (Geo == 3)
            #SPHERE
            x += ds*sqrt(1-mu^2)*cos(phi)
            y += ds*sqrt(1-mu^2)*sin(phi)
        end
    end
    return (x,y,z)
end

function getZone(x,y,z,low_edges,high_edges)
    r = sqrt(x^2 + y^2 + z^2)
    return argmax((r.>low_edges).*(r .<= high_edges))
end

function getNextZone(Geo,x,y,z,mu,phi,ds,low_edges,high_edges)
    x,y,z = updatePosition(Geo,x,y,z,mu,phi,ds)
    r = sqrt(x^2 + y^2 + z^2)
    return argmax((r.>low_edges).*(r .<= high_edges))
end

function getNextRadius(Geo,x,y,z,mu,phi,ds)
    x,y,z = updatePosition(Geo,x,y,z,mu,phi,ds)
    r = sqrt(x^2 + y^2 + z^2)
    return r
end

function getRadius(x,y,z)
    return sqrt(x^2 + y^2 + z^2)
end

function getPath(Geo,x,y,z,mu,phi,low_edges,high_edges,Nx)

    zone = getZone(x,y,z,low_edges,high_edges)

    if (Geo == 1)
        # direction to sweep
        if (mu >= 0)
            zoneList = (zone):1:Nx # going right
        else
            zoneList = (zone):-1:1 # going left
        end
        # partial zone
        ds = distance_to_edge(Geo,x,y,z,mu,phi,low_edges,high_edges,zone)
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
            ds = distance_to_edge(Geo,x,y,z,mu,phi,low_edges,high_edges,zone) + 1e-9
            append!(dsList, ds)
            append!(zoneList, zone)
            x,y,z = updatePosition(Geo,x,y,z,mu,phi,ds)
            if ((zone == Nx) && (ds == Inf))
                alive = false
            end
            if (getRadius(x,y,z) >= high_edges[Nx])
                alive = false
            end
        end
    end
    return (zoneList, dsList)
end

function getDim(Geo, hasLeft, hasRight)
    Dim = 0
    if (Geo == 1)
        Dim += 2 # volumetric source dimensions
        if (hasLeft)
            Dim += 2
        end
        if (hasRight)
            Dim += 2
        end
    else
        Dim += 3 # volumetric source dimensions
        if (hasLeft)
            Dim += 3
        end
        if (hasRight)
            Dim += 3
        end
    end
    return Dim
end

function inf_Med_BC_Mu(rn, x, siga)
    # parameters
    siga = siga[1]
    q0 = 0.75
    q1 = 0.5

    denominator = q0/(2*siga) + x*q1/(2*siga) - q1/(3*siga^2)
    f(mu) = (mu^2*q0/(2*siga) + mu^2*x*q1/(2*siga) - mu^3*q1/(3*siga^2))/denominator - rn

    a = 0
    b = 1
    N=0
    c = 0
    maxiter = 100
    tol = 1e-8

    # bisection method
    while (N<=maxiter) && ((b-a)*0.5 > tol)
        c = (a+b)*0.5
        if ((f(c)<0) && (f(a)<0)) || ((f(c)>0) && (f(a)>0))
            a = c
        else
            b = c
        end
        N += 1
    end
    return c
end
