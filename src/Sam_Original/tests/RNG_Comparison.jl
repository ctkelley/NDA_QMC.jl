"""
Sam Pasmann
"""
###############################################################################
#### Packages/Functions
###############################################################################

include("qmc_init.jl")
include("qmc_sweep.jl")
include("qmc_source_iteration.jl")
using PyPlot
using DelimitedFiles
using LinearAlgebra
using Statistics
pygui(true)

###############################################################################
#### Import Cross Section Data
###############################################################################

D = readdlm(joinpath(@__DIR__, "HDPE_data/D_12G_HDPE.csv"), ',', Float64)
siga = readdlm(joinpath(@__DIR__, "HDPE_data/Siga_12G_HDPE.csv"), ',', Float64)
sigs = readdlm(joinpath(@__DIR__, "HDPE_data/Scat_12G_HDPE.csv"), ',', Float64)
sigs = sigs
sigt = 1 ./ (3*D)
#sigt = siga + sum(sigs,dims=2)

###############################################################################
#### Parameters
###############################################################################

Nx = 50 #number of tally cells
na2 = 11 #number of angles for angular mesh
s = [1] #parameter in Garcia/Siewert
nTimes = 30
N = [2^3, 2^6, 2^9, 2^12]
sobol = zeros(size(N)[1], nTimes)
golden = zeros(size(N)[1], nTimes)
random = zeros(size(N)[1], nTimes)
#sigt = [1]
#sigs = [0.8]

###############################################################################
#### Source Iteration Call
###############################################################################

# analytical solution
source_strength = 1.0
G = size(sigt)[1] # number of groups
Q = source_strength*ones(G) # source
flux = inv(Diagonal(sigt[:,1]) .- sigs)*Q # returns diagonal matrix

for i in 1:size(N)[1]
    for j in 1:nTimes
        # function calls
        rng = "Golden"
        qmc_data = qmc_init(N[i], Nx, na2, s, sigs, sigt, rng)
        phi_avg_Golden, phi_edge, dphi, J_avg, J_edge, history, itt = qmc_source_iteration(s,qmc_data)

        rng = "Sobol"
        qmc_data = qmc_init(N[i], Nx, na2, s, sigs, sigt, rng)
        phi_avg_Sobol, phi_edge, dphi, J_avg, J_edge, history, itt = qmc_source_iteration(s,qmc_data)

        rng = "Random"
        qmc_data = qmc_init(N[i], Nx, na2, s, sigs, sigt, rng)
        phi_avg_Random, phi_edge, dphi, J_avg, J_edge, history, itt = qmc_source_iteration(s,qmc_data)

        sobol[i,j] = sum(abs.(mean(phi_avg_Sobol, dims=1)'[:,1]-flux))
        golden[i,j] = sum(abs.(mean(phi_avg_Golden, dims=1)'[:,1]-flux))
        random[i,j] = sum(abs.(mean(phi_avg_Random, dims=1)'[:,1]-flux))
    end
end


###############################################################################
#### Plots
###############################################################################
midpoints = qmc_data.midpoints

figure()
title("Mean Absolute Error")
plot(N, mean(sobol, dims=2), label="Sobol")
plot(N, mean(golden, dims=2), label="Golden")
plot(N, mean(random, dims=2), label="Random")
xscale("log")
yscale("log")
ylabel("Mean Absolute Error")
xlabel("Number of Particles")
legend()
