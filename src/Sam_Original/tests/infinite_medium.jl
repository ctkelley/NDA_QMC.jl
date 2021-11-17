###############################################################################
#### Packages/Functions
###############################################################################

include("../qmc_init.jl")
include("../qmc_sweep.jl")
include("../qmc_source_iteration.jl")
using PyPlot
pygui(true)

###############################################################################
#### Parameters
###############################################################################

Nx = 40     # number of tally cells
na2 = 11    # number of angles for angular mesh
s = [1]     # parameter in Garcia/Siewert
N = 2^13    # number of particles per source itertion
LB = 0      # left bound
RB = 2      # right bound
sigt = [2.0]
sigs = 0.5*sigt
siga = sigt - sigs
geometry = "Slab"   # infinite medium problems are only configured for slab geometries right now
generator = "Sobol"

###############################################################################
#### Function Call
###############################################################################

qmc_data = linear_infMed_init(geometry, generator, N, LB, RB, Nx, na2, s, sigs, sigt)
phi_avg, phi_edge, dphi, J_avg, J_edge, history, itt = qmc_source_iteration(s,qmc_data)

###############################################################################
#### Plotting
###############################################################################
midpoints = qmc_data.midpoints
flux = qmc_data.phi_avg
sol = qmc_data.phi_true

figure()
plot(midpoints, flux, label = generator)
plot(midpoints, sol, label= "solution")
ylabel("Cell Average Flux")
xlabel("Cell Midpoints")
title("Scalar Flux")
legend()
