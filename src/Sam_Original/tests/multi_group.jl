"""
Sam Pasmann
"""
###############################################################################
#### Packages/Functions
###############################################################################

include("../qmc_init.jl")
include("../qmc_sweep.jl")
include("../qmc_source_iteration.jl")
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Statistics
pygui(true)

###############################################################################
#### Import Cross Section Data
###############################################################################

D = readdlm(joinpath(@__DIR__, "../HDPE_data/D_12G_HDPE.csv"), ',', Float64)
siga = readdlm(joinpath(@__DIR__, "../HDPE_data/Siga_12G_HDPE.csv"), ',', Float64)
sigs = readdlm(joinpath(@__DIR__, "../HDPE_data/Scat_12G_HDPE.csv"), ',', Float64)
sigs = reverse(sigs, dims=2)
sigt = 1 ./ (3*D)

###############################################################################
#### Parameters
###############################################################################

Nx = 50     # number of tally cells
na2 = 11    # number of angles for angular mesh
s = [1]     # parameter in Garcia/Siewert
N = 2^11   # number of particles per itertion per source
LB = 0      # left bound
RB = 500      # right bound
geometry = "Slab"
generator = "Sobol"

qmc_data = qmc_init(geometry, generator, N, LB, RB, Nx, na2, s, sigs, sigt)
phi_avg_Sobol, phi_edge, dphi, J_avg, J_edge, history, itt = qmc_source_iteration(s,qmc_data)

###############################################################################
#### Analytic Solution
###############################################################################
source_strength = 1.0
G = size(sigt)[1] # number of groups
Q = source_strength*ones(G) # source
flux = inv(Diagonal(sigt[:,1]) .- sigs)*Q # returns diagonal matrix

###############################################################################
#### Plots
###############################################################################
midpoints = qmc_data.midpoints

figure()
title("Scalar Flux")
plot(midpoints, sum(phi_avg_Sobol, dims=2), label="Sobol")
#plot(midpoints, sum(phi_avg_Golden, dims=2), label="Golden")
#plot(midpoints, sum(phi_avg_Random, dims=2), label="Random")
#plot(midpoints, sum(flux)*ones(Nx), label="Analytic Sol")
ylabel("cell averaged flux")
xlabel("midpoints")
legend()


counter = 1
figure()
color_sequence = ["#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
                  "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
                  "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
                  "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"]
title("Group Scalar Flux")
for i in 1:1:G
    plot(midpoints, phi_avg_Sobol[:,i], label=i,color=color_sequence[i] ) #color=color_sequence[i]
    plot(midpoints,flux[i]*ones(Nx),"--",color=color_sequence[i] )
end
#plot(midpoints,phi_avg, label="QMC")
ylabel("cell averaged flux")
xlabel("midpoints")
legend()