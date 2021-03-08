include("QMC_sweep.jl")

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

#cell averaged scalar flux
figure(figsize = (5,10))

subplot(311)
plot(midpoints,phi_avg)
ylabel("phi")
xlabel("midpoints")

subplot(312)
plot(midpoints,dphi)
ylabel("dphi")
xlabel("midpoints")

subplot(313)
plot(midpoints,J_avg)
ylabel("current")
xlabel("midpoints")
#yscale("log")

#left exit bins
figure(figsize = (5,10))
subplot(211)
gsSol = zeros((11,2))
gsSol[:,1] = -1*[0.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]
gsSol[:,2] = [0.58966,0.53112,0.44328,0.38031,0.33296,0.29609,0.26656,0.24239,
                                                        0.22223,0.20517,0.19055]
plot(gsSol[:,1], gsSol[:,2],"+")
plot(exit_left_bins[:,1], exit_left_bins[:,2])

#right exit bins
subplot(212)
gsSol = zeros((11,2))
gsSol[:,1] = 1*[0.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]
gsSol[:,2] = [6.08E-06,6.93E-06,9.64E-06,1.62E-05,4.39E-05,1.69E-04,5.73E-04,
                                            1.51E-03,3.24E-03,5.96E-03,9.77E-03]
plot(gsSol[:,1], gsSol[:,2],"+")
plot(exit_right_bins[:,1], exit_right_bins[:,2])
