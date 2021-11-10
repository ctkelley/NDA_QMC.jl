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
N = 2^11    # number of particles per source itertion
LB = 0      # left bound
RB = 10      # right bound
#sigt = [6.161861606013174]
#sigs = [5.719904403604628]
sigt = [1.0]
sigs = 0.5*sigt
siga = sigt - sigs
geometry = "Slab"
generator = "Sobol"

###############################################################################
#### Function Call
###############################################################################

qmc_data = infMed_init(geometry, generator, N, LB, RB, Nx, na2, s, sigs, sigt)
@time begin
phi_avg, phi_edge, dphi, J_avg, J_edge, history, itt = qmc_source_iteration(s,qmc_data)
end

###############################################################################
#### Plotting
###############################################################################
midpoints = qmc_data.midpoints
flux = qmc_data.phi_avg

figure()
plot(midpoints, flux, label = generator)
ylabel("Cell Average Flux")
xlabel("Cell Midpoints")
title("Scalar Flux")
legend()
