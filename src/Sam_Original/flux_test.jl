"""
Sam Pasmann
07/14/2021

Scalar Flux at cell edges appear to have a bug in a cell somewhere near the
center which has a large spike with increasing N.

Update:
    bug only showing up when N=10^4

    The Sobol sequence casts nets of 2^n . so N should be a power of 2

    For some reason even if N is a power of 2, skipping can still cause spikes
"""


include("qmc_init.jl")
include("qmc_sweep.jl")
using PyPlot
#using DelimitedFiles
using CSV
pygui(true)

###############################################################################
#### Import Cross Section Data
###############################################################################

D = readdlm(joinpath(@__DIR__, "HDPE_data/D_12G_HDPE.csv"), ',', Float64)
siga = readdlm(joinpath(@__DIR__, "HDPE_data/Siga_12G_HDPE.csv"), ',', Float64)
sigs = readdlm(joinpath(@__DIR__, "HDPE_data/Scat_12G_HDPE.csv"), ',', Float64)
#sigt = 1 ./ (3*D)
sigt = siga + sum(sigs,dims=2)

###############################################################################
#### Parameters
###############################################################################

N = [10^2, 10^3, 10^4, 10^5] #number of QMC particles
N = [2^8, 2^11, 2^14, 2^17]
Nx = 100 #number of tally cells
na2 = 11 #number of angles for angular mesh
s = [1] #parameter in Garcia/Siewert
figure(1, figsize = (20,4))
###############################################################################
#### Function Call
###############################################################################
for i in 1:length(N)
    # initialize qmc
    qmc_data = qmc_init(N[i], Nx, na2, s, sigs, sigt)
    # call qmc SI
    phi_avg, phi_edge, dphi, J_avg, J_edge = qmc_sweep(qmc_data)

    # Plotting
    midpoints = qmc_data.midpoints
    edges = qmc_data.edges
    #edge flux
    num = 150 + i
    subplot(num)
    suptitle("Average Scalar Flux")
    plot(midpoints,phi_avg)
    ylabel("cell averaged flux")
    xlabel("midpoints")
    title("N = $(N[i])")
end
