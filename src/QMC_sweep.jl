using PyPlot
using Random
using SpecialFunctions
using Sobol
using LinearAlgebra
using JLD2
pygui(true)

###############################################################################
# Documentation
###############################################################################
"""
1D, 1 group, slab, Quasi Monte Carlo functions. Particularly designed for the
problems outlined in Garcia et al. where Sigma_s = exp(x/s) and the angular
flux distribution is tallied at the boundaries.

move_part:
        Moves particle from initial position to time or slab boundary, tallying
        scalar flux and current per cell.
QMC_Sweep:
        Initializes variables, calls move_part for each particle.


    Inputs:
    -------
        source_strength:
            Scalar value, strength of internal source
        N:
            Number of internal source particles
        has_right:
            True or false, source right boundary
        has_left:
            True or false, source left boundary
        Nb:
            Number of boundary source particles
        Nx:
            Number of scalar flux and current tallying cells
        Lx:
            Width of slab (cm)
        nbins:
            Number of exit angular flux tally bins
        dt:
            Time increment (s). Set dt = Inf if S.S. is desired
        v:
            Particle velocity (cm/s)
        sigt:
            Scalar value, total macroscopic cross section
        c:
            Scattering cross section varies with exp(x/c)
        phi_avg:
            Optional input. Defaults to zeros(Nx)


    Outputs:
    --------
        midpoints:
            Array of length Nx, midpoints of cells
        edges:
            Array of length Nx+1, cell boundaries
        phi_avg:
            Array of length Nx, average scalar flux across each cell
        phi_edge:
            Array of length Nx+1, scalar flux at cell boundaries
        dphi:
            Spatial derivative of phi_avg
        J_avg:
            Array of length Nx, average cell current
        J_edge:
            Array of length Nx+1, current at cell boundaries
        exit_right_bins:
            Array of length nbins, angular flux distribution at right boundary
        exit_left_bins:
            Array of length nbins, angular flux distribution at left boundary
"""

###############################################################################
# Functions
###############################################################################

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
            phi_avg[zone] += weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone] #implicit capture for track-length
            phi_s[zone] += (weight*exp(-x/c))*(1 - exp(-(ds*(sigt[zone] + mu/c))))/(sigt[zone] + mu/c)/dxs[zone] #weighted average version
            J_avg[zone] += phi_avg[zone]*mu
            #J_avg[zone] += weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone]*mu
        else
            phi_avg[zone] += weight*ds/dxs[zone] #implicit capture for track-length
            phi_s[zone] += phi_avg[zone]
            J_avg[zone] += phi_avg[zone] *mu
        end
        weight *= exp(-(ds*sigt[zone])) #decrease the weight
        J_edge[zone+1] += weight
        phi_edge[zone+1] += weight/mu #current density tally
        s += ds

        # Full zones
        z_prop = zone+1
        while (s+ds_zone <= s_max) && (z_prop <= Nx)
            if (sigt[z_prop] > 1e-12)
                ds += ds_zone
                phi_avg[z_prop] += weight*(1 - exp(-(ds_zone*sigt[z_prop])))/sigt[z_prop]/dxs[z_prop] #implicit capture for track-length
                phi_s[z_prop] += (weight*exp(-low_edges[z_prop]/c))*(1 - exp(-(ds_zone*(sigt[z_prop] + mu/c))))/(sigt[z_prop] + mu/c)/dxs[z_prop]
                J_avg[z_prop] += phi_avg[z_prop]*mu
                #J_avg[z_prop] += weight*(1 - exp(-(ds*sigt[z_prop])))/sigt[z_prop]/dxs[z_prop]*mu
                s += ds_zone
            else
                phi_avg[z_prop] += weight*ds_zone/dxs[z_prop] #implicit capture for track-length
                phi_s[z_prop] += phi_avg[z_prop]
                J_avg[z_prop] += phi_avg[z_prop]*mu
                s += ds_zone
            end
            weight *= exp(-(ds_zone*sigt[z_prop]))
            J_edge[z_prop+1] += weight
            phi_edge[z_prop+1] += weight/mu
            z_prop += 1
        end

        # Last partial zone
        if (z_prop < Nx) && (s < s_max)
            zone = argmax(1*((s*mu).>=low_edges).*((s*mu) .< high_edges))
            ds = s_max-s
            if (sigt[zone] > 1e-12)
                phi_avg[zone] += weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone] #implicit capture for track-length
                phi_s[zone] += (weight*exp(-x/c))*(1 - exp(-(ds*(sigt[zone] + mu/c))))/(sigt[zone] + mu/c)/dxs[zone] #weighted average version
                J_avg[zone] += phi_avg[zone]*mu
                #J_avg[zone] += weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone]*mu
            else
                phi_avg[zone] += weight*ds/dxs[zone] #implicit capture for track-length
                phi_s[zone] += phi_avg[zone]
                J_avg[zone] += phi_avg[zone]*abs(mu)
            end
            weight *= exp(-(ds*sigt[zone])) #decrease the weight
            J_edge[zone+1] += weight #current density tally
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
            phi_avg[zone] += weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone] #implicit capture for track-length
            phi_s[zone] += (weight*exp(-x/c))*(1 - exp(-(ds*(sigt[zone] + abs(mu)/c))))/(sigt[zone] + abs(mu)/c)/dxs[zone]
            J_avg[zone] += phi_avg[zone]*abs(mu)
            #J_avg[zone] += weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone]*abs(mu)
        else
            phi_avg[zone] += weight*ds/dxs[zone]
            phi_s[zone] += phi_avg[zone]
            J_avg[zone] += phi_avg[zone]*abs(mu)
        end
        weight *= exp.(-(ds*sigt[zone]))
        J_edge[zone] -= weight
        phi_edge[zone] -= weight/abs(mu)
        s += ds
        # Full zones
        z_prop = zone-1
        while (s+ds_zone <= s_max) && (z_prop >= 1)
            if (sigt[z_prop] > 1e-12)
                phi_avg[z_prop] += weight*(1 - exp(-(ds_zone*sigt[z_prop])))/sigt[z_prop]/dxs[z_prop] #implicit capture for track-length
                phi_s[z_prop] += (weight*exp(-high_edges[z_prop]/c))*(1 - exp(-(ds_zone*(sigt[z_prop] + abs(mu)/c))))/(sigt[z_prop] + abs(mu)/c)/dxs[z_prop]
                J_avg[z_prop] += phi_avg[z_prop]*abs(mu)
                #J_avg[z_prop] += weight*(1 - exp(-(ds*sigt[z_prop])))/sigt[z_prop]/dxs[z_prop]*abs(mu)
                s += ds_zone
            else
                phi_avg[z_prop] += weight*ds_zone/dxs[z_prop] #implicit capture for track-length
                phi_s[z_prop] += phi_avg[z_prop]
                J_avg[z_prop] += phi_avg[z_prop]*abs(mu)
                s += ds_zone
            end
           weight *= exp(-(ds_zone*(sigt[z_prop])))
           J_edge[z_prop] -= weight
           phi_edge[z_prop] -= weight/abs(mu)
           z_prop -= 1
       end

        # Last partial zone
        if (z_prop > 1) && (s < s_max)
            zone = argmax(1*((s*mu).>=low_edges).*((s*mu) .< high_edges))
            ds = (low_edges[zone] - s*mu)/mu
            if (sigt[zone] > 1e-12)
                phi_avg[zone] += weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone] #implicit capture for track-length
                phi_s[zone] += (weight*exp(-x/c))*(1 - exp(-(ds*(sigt[zone] + abs(mu)/c))))/(sigt[zone] + abs(mu)/c)/dxs[zone]
                J_avg[zone] += phi_avg[zone]*abs(mu)
                #J_avg[zone] += weight*(1 - exp(-(ds*sigt[zone])))/sigt[zone]/dxs[zone]*abs(mu)
            else
                phi_edge[zone] += weight*ds/dxs[zone] #implicit capture for track-length
                phi_s[zone] += phi_edge[zone]
                J_avg[zone] += phi_edge[zone]*abs(mu)
            end
            weight *= exp.(-(ds*sigt[zone]))
            J_edge[zone] -= weight
            phi_edge[zone] -= weight/abs(mu)
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


function QMC_Sweep(N,Nb,Nx,Lx,nbins,dt,v,sigt,source,has_left,has_right,c, phi_avg = zeros(Nx))

    rng = SobolSeq(3)
    if (N >0)
        skip(rng,N) #skipping the expected number is suggested for Sobol
    end

    rng_bndl = 0
    if (has_left > 0)
        rng_bndl = SobolSeq(2)
        skip(rng_bndl,Nb) #skipping the expected number is suggested for Sobol
    end
    rng_bndr = 0
    if (has_right > 0)
        rng_bndr = SobolSeq(2)
        skip(rng_bndr,Nb) #skipping the expected number is suggested for Sobol
    end

    #define tally mesh
    dx = Lx/Nx
    low_edges = range(0, stop=Lx-dx, length=Nx)
    high_edges = low_edges.+dx
    dxs = high_edges - low_edges
    midpoints = 0.5*(high_edges + low_edges)
    edges = range(0, stop=Lx, length=Nx+1)

    fl(mu) =  sqrt(mu)
    fr(mu) = -sqrt(mu)

    #define angular flux mesh
    #exit_left_bins data structure to hold the exiting angular flux,
    #the first column has the bin centers and the
    #and the second holds the values.
    #exit_right_bins is the same
    dmu = 1/nbins
    #right bins
    exit_right_bins = zeros((nbins,2))
    exit_right_bins[:,1] = range(dmu/2, stop = 1- dmu/2, step = dmu)
    exit_right_bins[:,2] .= 0
    #left bins
    exit_left_bins = zeros((nbins,2))
    exit_left_bins[:,1] = -1*range(dmu/2, stop = 1- dmu/2, step = dmu)
    exit_left_bins[:,2] .= 0

    #phi_avg is defaulted to = zeros(Nx)
    phi_edge = zeros(Nx+1)
    dphi = zeros(Nx)
    phi_s = zeros(Nx) .+ 1e-6

    J_avg = zeros(Nx)
    J_edge = zeros(Nx+1)

    sigt = ones(Nx)*sigt
    source = source_strength*ones(Nx)
    sigsFunc(x) = exp.(-x/c)
    sigs = sigsFunc(midpoints)
    q = phi_avg.*sigs .+ source

    #do boundary sources
    if (has_left >0)
        for i in 1:Nb
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
    end
    if (has_right >0)
        for i in 1:Nb
            tmp_rnd = next!(rng_bndr) #for Sobol
            x = Lx #start at right boundary
            mu = fr(tmp_rnd[1]) #mu in -1 to 0
            if dt == Inf
                t = 0
            else
                t = dt*fr(tmp_rnd[2]) # start time somewhere in dt
            end
            #determine zone (right boundary)
            zone = Nx
            #compute weight
            weight = 1/Nb*0.5
            #set phi_edge
            phi_edge[Nx] = 1
            #how far do I travel when I cross a zone
            ds_zone = dx/abs(mu)
            J_edge[zone] -= weight
            move_part(mu,zone,x,Nx,dt,t,v,high_edges,low_edges,dxs,weight,ds_zone,
                                    phi_avg, phi_edge, phi_s, J_avg, J_edge,sigt,
                                            exit_right_bins,exit_left_bins,c)
        end
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

   return (midpoints,edges, phi_avg, phi_edge, dphi, J_avg, J_edge,
                                                exit_right_bins, exit_left_bins)
end
