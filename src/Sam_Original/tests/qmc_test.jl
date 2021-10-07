using PyPlot
include("qmc_init.jl")
include("qmc_source_iteration.jl")

###############################################################################
#### Inputs
###############################################################################
N = 10^4 #number of QMC particles
Nx = 50 #number of tally cells
na2 = 11 #number of angles for angular mesh
s = 1 #parameter in Garcia/Siewert

#phi_edge = qmc_data.phi_edge
###############################################################################
#### Function Call
###############################################################################
# initialize qmc
qmc_data = qmc_init(N, Nx, na2, s)
# call qmc SI
phi_avg, phi_edge, dphi, J_avg, J_edge, psi_right, psi_left, history, itt = qmc_source_iteration(s,qmc_data)


###############################################################################
#### Plotting
###############################################################################
####
iteration = 1:itt

midpoints = qmc_data.midpoints
edges = qmc_data.edges

####
figure(figsize = (6,12))

subplot(311)
plot(midpoints,phi_avg)
ylabel("phi midpoints")
xlabel("midpoints")
title("Cell Averaged Scalar Flux")

# cell Averaged Spatial Derivative of Scalar Flux
subplot(312)
plot(midpoints,dphi)
ylabel("dphi")
xlabel("cell midpoints")
title("Cell Flux Derivative")

# cell averaged current
subplot(313)
plot(midpoints,J_avg)
ylabel("J avg")
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
plot(psi_left[:,1], psi_left[:,2],label="QMC_sweep")
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
plot(psi_right[:,1], psi_right[:,2], label="QMC_sweep")
xlabel("mu")
ylabel("Right Boundary - Angular Flux Dist.")
legend()
