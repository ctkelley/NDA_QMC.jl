"""
    rngInit(generator, N, dim)
Random Number Generator (rng) Initilizer.
...
# Arguments
* `generator::string`: specify rng type. "Random", "Sobol", "Golden".
* `N::Integer`: number of particles.
* `dim::Integer`: total number of dimensions.
...
...
# Outputs
* `rng::Array{Float64,2}`: NxDim matrix of random numbers from specified rng.
...
"""
function rngInit(generator, N, dim)
    #   skipping the expected number is suggested for Sobol
    #   but has been causing spikes for higher particle counts

    # need to  create a matrix of random numbers by number of dimensions
    # the sobol package makes this difficult
    if (generator == "Sobol")
        sobol = SobolSeq(dim)
        rng = zeros(N, dim)
        for i in 1:N
            temp = next!(sobol)
            rng[i,:] = temp
        end
    elseif (generator == "Random")
        rng = rand(N,dim)
    elseif (generator == "Golden")
        rng = GoldenSequence(dim)
    else
        error("RNG must be 'Sobol', 'Golden', or 'Random'. ")
    end
    return rng
end

"""
    nextBoundaryRN(rng, i, generator, Geo, Dim)
Take next two random numbers from rng required for boundary source.
randPhi will return as zero if in slab geometry.

...
# Arguments
* `rng::Array{Float64,2}`: NxDim matrix of random numbers.
* `i::Float64`: Particle number.
* `generator::string`: specify rng type. "Random", "Sobol", "Golden".
* `Geo::Integer`: geometry specification.
* `Dim::Integer`: total number of dimensions.
...
...
# Outputs
* `randMu::Float64`: random number between [0,1] for Mu.
* `randPhi::Float64`: random number between [0,1] for Phi. Zeros for slab geometry.
...
"""
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

"""
    nextRN(rng, i, generator, Geo, Dim)
Take next three random numbers from rng required for volumetric source.
randPhi will return as zero if in slab geometry.

...
# Arguments
* `rng::Array{Float64,2}`: NxDim matrix of random numbers.
* `i::Float64`: Particle number.
* `generator::string`: specify rng type. "Random", "Sobol", "Golden".
* `Geo::Integer`: geometry specification.
* `Dim::Integer`: total number of dimensions.
...
...
# Outputs
* `randX::Float64`: random number between [0,1] for X postion.
* `randMu::Float64`: random number between [0,1] for Mu angle.
* `randPhi::Float64`: random number between [0,1] for Phi angle. Zero for slab geometry.
...
"""
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

"""
    cellVolume(Geo, zone, low_edges, high_edges)
Calculate cell volume given mesh and geometry type.

...
# Arguments
* `Geo::Integer`: geometry specification.
* `zone::Integer`: current zone of particle.
* `low_edges::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`: lower boundaries of cells.
* `high_edges::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`: upper boundaries of cells.
...
...
# Outputs
* `Float64`: volume of current cell.
...
"""
function cellVolume(Geo, zone, low_edges, high_edges)
    if (Geo == 1) # slab
        return (high_edges[zone] - low_edges[zone])
    elseif (Geo == 2) # cylinder
        return pi*(high_edges[zone]^2 - low_edges[zone]^2)
    elseif (Geo == 3) # sphere
        return (4/3)*pi*(high_edges[zone]^3 - low_edges[zone]^3)
    end
end

"""
    surfaceArea(Geo, r)
Surface area by geometry for boundary source adjustment.

...
# Arguments
* `Geo::Integer`: geometry specification.
* `r::Float64`: radius.
...
...
# Outputs
* `Float64`: surface area given geometry type.
...
"""
function surfaceArea(Geo, r)
    if (Geo == 1) # slab
        return (0.5)
    elseif (Geo == 2) # cylinder
        return 2*pi*(r)
    elseif (Geo == 3) # sphere
        return (4)*pi*(r^2)
    end
end

"""
    sphere_edge(x,y,z,mu,phi,edge)
Solve spherical surface interation quadratic.

...
# Arguments
* `x::Float64`: x position.
* `y::Float64`: y position.
* `z::Float64`: z position.
* `mu::Float64`: mu angle.
* `phi::Float64`: phi angle.
* `edge::Float64`: radius of surface.
...
...
# Outputs
* `ds::Float64`: distance to edge. Will return `Inf` if no interation occurs.
...
"""
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

"""
    cylinder_edge(x,y,mu,phi,edge)
Solve cylindrical surface interation quadratic.

...
# Arguments
* `x::Float64`: x position.
* `y::Float64`: y position.
* `mu::Float64`: mu angle.
* `phi::Float64`: phi angle.
* `edge::Float64`: radius of surface.
...
...
# Outputs
* `ds::Float64`: distance to edge. Will return `Inf` if no interation occurs.
...
"""
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

"""
    slab_edge(mu,x,edge)
Solve distance to edge in a slab.

...
# Arguments
* `mu::Float64`: mu angle.
* `x::Float64`: x position.
* `edge::Float64`: radius of surface.
...
...
# Outputs
* `Float64`: distance to edge.
...
"""
function slab_edge(mu,x,edge)
    return (edge - x)/mu
end

"""
    distance_to_edge(Geo,x,y,z,mu,phi,low_edges,high_edges,zone)
Determine which of the two surfaces the particle is traveling to and
return distance to that surface.

...
# Arguments
* `Geo::Integer`: geometry type.
* `x::Float64`: x position.
* `y::Float64`: y position.
* `z::Float64`: z position.
* `mu::Float64`: mu angle.
* `phi::Float64`: phi angle.
* `low_edges::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`: lower boundaries of cells.
* `high_edges::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`: upper boundaries of cells.
* `zone::Integer`: current cell.
...
...
# Outputs
* `ds::Float64`: distance to edge.
...
"""
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
    return (ds)
end

"""
    updatePosition(Geo,x,y,z,mu,phi,ds)
Update x, y, and z positions given angles and distance traveled (ds).

...
# Arguments
* `Geo::Integer`: geometry type.
* `x::Float64`: x position.
* `y::Float64`: y position.
* `z::Float64`: z position.
* `mu::Float64`: mu angle.
* `phi::Float64`: phi angle.
* `ds::Float64`: distance traveled.
...
...
# Outputs
* `x::Float64`: new x position.
* `y::Float64`: new y position.
* `z::Float64`: new z position.
...
"""
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

"""
    getZone(x,y,z,low_edges,high_edges)
Return current zone of partilce.

...
# Arguments
* `x::Float64`: x position.
* `y::Float64`: y position.
* `z::Float64`: z position.
* `low_edges::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`: lower boundaries of cells.
* `high_edges::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`: upper boundaries of cells.
...
...
# Outputs
* `zone::Integer`: current cell.
...
"""
function getZone(x,y,z,low_edges,high_edges)
    r = sqrt(x^2 + y^2 + z^2)
    return argmax((r.>low_edges).*(r .<= high_edges))
end

"""
    getNextZone(Geo,x,y,z,mu,phi,ds,low_edges,high_edges)
Return the next zone the particle will be in given current angle and distance.

...
# Arguments
* `Geo::Integer`: geometry type.
* `x::Float64`: x position.
* `y::Float64`: y position.
* `z::Float64`: z position.
* `mu::Float64`: mu angle.
* `phi::Float64`: phi angle.
* `ds::Float64`: path length.
* `low_edges::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`: lower boundaries of cells.
* `high_edges::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`: upper boundaries of cells.
...
...
# Outputs
* `zone::Integer`: next cell.
...
"""
function getNextZone(Geo,x,y,z,mu,phi,ds,low_edges,high_edges)
    x,y,z = updatePosition(Geo,x,y,z,mu,phi,ds)
    r = sqrt(x^2 + y^2 + z^2)
    return argmax((r.>low_edges).*(r .<= high_edges))
end

"""
    getNextRadius(Geo,x,y,z,mu,phi,ds)
Return the next radius of the particle given current angle and distance.
...
# Arguments
* `Geo::Integer`: geometry type.
* `x::Float64`: x position.
* `y::Float64`: y position.
* `z::Float64`: z position.
* `mu::Float64`: mu angle.
* `phi::Float64`: phi angle.
* `ds::Float64`: path length.
...
...
# Outputs
* `r::Float64`: next radius.
...
"""
function getNextRadius(Geo,x,y,z,mu,phi,ds)
    x,y,z = updatePosition(Geo,x,y,z,mu,phi,ds)
    r = sqrt(x^2 + y^2 + z^2)
    return r
end

"""
    getRadius(x,y,z)
Return current radius of particle
...
# Arguments
* `x::Float64`: x position.
* `y::Float64`: y position.
* `z::Float64`: z position.
...
...
# Outputs
* `r::Float64`: radius.
...
"""
function getRadius(x,y,z)
    return sqrt(x^2 + y^2 + z^2)
end

"""
    getPath(Geo,x,y,z,mu,phi,low_edges,high_edges,Nx)
Generate a list of zones and path length(s) per zone that the particle will travel
from emission to exiting the volume.

In slab geometry this is a straigtfoward calculation, the particle will travel
through all zones in one direction (left or right) and will have the same path
length in each zone (minus the emission zone which will be a partial length).

In cylindrical and spherical geometries, the only way to determine which boundary
the particle will hit next (inner or outter radius) is to solve the quadratic for
both surfaces and take the shortest if both exist (one will always exist).

...
# Arguments
* `Geo::Integer`: geometry type.
* `x::Float64`: x position.
* `y::Float64`: y position.
* `z::Float64`: z position.
* `mu::Float64`: mu angle.
* `phi::Float64`: phi angle.
* `ds::Float64`: path length.
* `low_edges::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`: lower boundaries of cells.
* `high_edges::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`: upper boundaries of cells.
...
...
# Outputs
* `zoneList::Array{Float64,2}`:
* `dsList::Array{Float64,2}`:
...
"""
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

"""
    getDim(Geo, hasLeft, hasRight)
The dimensionality of the problem is needed to initialize the random number matrix.
The total number of dimensions is determined by the geometry and the presence of boundary
sources.

...
# Arguments
* `Geo::Integer`: geometry type.
* `hasLeft::Bool`: left boundary source present.
* `hasRight::Bool`: right boundary source present.
...
...
# Outputs
* `Dim::Integer`: total number of dimensions.
...
"""
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



"""
    inf_Med_BC_Mu(rn,x,siga)

Linear anisotropic direct inversion sampling of source term for angle.

Function is currently only built for one-group problems. Quadratic
source terms will require a separate function solve.

...
# Arguments
* `rn::Float64`: random number [0,1].
* `x::Float64`: spatial position.
* `siga::Array{Float64,1}`: absorption cross section
...
...
# Outputs
* `mu::Float64`: anisotropic sampled anlge.
...
"""
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
    mu = 0
    maxiter = 100
    tol = 1e-8

    # bisection method
    while (N<=maxiter) && ((b-a)*0.5 > tol)
        mu = (a+b)*0.5
        if ((f(mu)<0) && (f(a)<0)) || ((f(mu)>0) && (f(a)>0))
            a = mu
        else
            b = mu
        end
        N += 1
    end
    return mu
end
