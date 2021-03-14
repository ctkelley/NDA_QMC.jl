include("QMC_sweep.jl")

###############################################################################
# Function Documentation
###############################################################################
"""
QMC_Sweep:
        Performs one QMC sweep of a 1D, 1 group, slab, problem. Particularly 
        designed for the problems outlined in Garcia et al. where 
        Sigma_s = exp(-x/c) and the angular flux distribution is tallied at the 
        boundaries.

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
            Scattering cross section varies with exp(-x/c)
        phi_avg:
            Optional input of cell averaged scalar flux.
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
# Example Inputs
###############################################################################

dt = Inf# Inf #time increment
v = 1 #velocity (cm/s)
c = 1 #scattering exp(x/c) term from garica problem
Nx = 50 #number of tally bins
Lx = 5 #length of slab (cm)
nbins = 10 #number of exit angular flux tally bins

# left or right boundary source
has_left = true
has_right = false
#particles
Nb = 10^4 #boundary
N = 10^4 #internal
source_strength = 0.0 #internal source
sigt = 1 #total macroscopic cross section

#call the function
midpoints,edges, phi_avg,phi_edge,dphi,J_avg,J_edge, exit_right_bins, exit_left_bins= QMC_Sweep(N,Nb,Nx,Lx, nbins,dt,v,sigt,source_strength,has_left,has_right, c)

###############################################################################
# Example Plotting
###############################################################################
"""
Figure 1: 
          - cell averaged scalar flux 
          - cell averaged spatial derivative of scalar flux
          - cell averaged current 

Figure 2: 
          - Left boundary angular flux distribution
          - right boundary angular flux distribution

"""

#cell averaged scalar flux
figure(figsize = (5,10))

subplot(311)
plot(midpoints,phi_avg)
ylabel("phi")
xlabel("midpoints")
title("Cell Averaged Scalar Flux")

# cell Averaged Spatial Derivative of Scalar Flux
subplot(312)
plot(midpoints,dphi)
ylabel("dphi")
xlabel("midpoints")
title("Cell Averaged Spatial Derivative of Scalar Flux")

# cell averaged current
subplot(313)
plot(midpoints,J_avg)
ylabel("current")
xlabel("midpoints")
title("Cell Averaged Current")

#left exit bins
figure(figsize = (5,10))
subplot(211)
gsSol = zeros((11,2))
gsSol[:,1] = -1*[0.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]
gsSol[:,2] = [0.58966,0.53112,0.44328,0.38031,0.33296,0.29609,0.26656,0.24239,
                                                        0.22223,0.20517,0.19055]
plot(gsSol[:,1], gsSol[:,2],"+",label="Garcia et al.")
plot(exit_left_bins[:,1], exit_left_bins[:,2],label="QMC_sweep")
xlabel("-mu")
ylabel("Left Boundary - Angular Flux Dist.")
title("Angular Flux Exit Distributions from Garcia et al., c = 1")
legend()

#right exit bins
subplot(212)
gsSol = zeros((11,2))
gsSol[:,1] = 1*[0.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]
gsSol[:,2] = [6.08E-06,6.93E-06,9.64E-06,1.62E-05,4.39E-05,1.69E-04,5.73E-04,
                                            1.51E-03,3.24E-03,5.96E-03,9.77E-03]
plot(gsSol[:,1], gsSol[:,2],"+", label="Garcia et al.")
plot(exit_right_bins[:,1], exit_right_bins[:,2], label="QMC_sweep")
xlabel("mu")
ylabel("Right Boundary - Angular Flux Dist.")
legend()
