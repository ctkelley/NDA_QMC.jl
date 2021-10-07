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
#sigt = siga + sum(sigs,dims=2)

###############################################################################
#### Parameters
###############################################################################

Nx = 50 #number of tally cells
na2 = 11 #number of angles for angular mesh
s = [1] #parameter in Garcia/Siewert
N = 2^11

###############################################################################
#### Source Iteration Call
###############################################################################

# function calls
rng = "Golden"
qmc_data = qmc_init(N, Nx, na2, s, sigs, sigt, rng)
phi_avg_Golden, phi_edge, dphi, J_avg, J_edge, history, itt = qmc_source_iteration(s,qmc_data)

rng = "Sobol"
qmc_data = qmc_init(N, Nx, na2, s, sigs, sigt, rng)
phi_avg_Sobol, phi_edge, dphi, J_avg, J_edge, history, itt = qmc_source_iteration(s,qmc_data)

rng = "Random"
qmc_data = qmc_init(N, Nx, na2, s, sigs, sigt, rng)
phi_avg_Random, phi_edge, dphi, J_avg, J_edge, history, itt = qmc_source_iteration(s,qmc_data)

# analytical solution
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
plot(midpoints, sum(phi_avg_Golden, dims=2), label="Golden")
plot(midpoints, sum(phi_avg_Random, dims=2), label="Random")
plot(midpoints, sum(flux)*ones(Nx), label="Analytic Sol")
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
    plot(midpoints, phi_avg_Sobol[:,i], label=i,color=color_sequence[counter] ) #color=color_sequence[i]
    plot(midpoints,flux[i]*ones(Nx),"--",color=color_sequence[counter] )
    counter += 1
end
#plot(midpoints,phi_avg, label="QMC")
ylabel("cell averaged flux")
xlabel("midpoints")
legend()



"""
figure()
plot(1:G, abs.(mean(phi_avg_Sobol, dims=1)'[:,1]-flux), label="Sobol")
#plot(1:G, abs.(mean(phi_avg_Golden, dims=1)'[:,1]-flux), label="Golden")
#plot(1:G, abs.(mean(phi_avg_Random, dims=1)'[:,1]-flux), label="Random")
ylabel("Mean Absolute Error")
xlabel("Group")
title("Average Group Error")
legend()
"""
###############################################################################
#### Single Function Call
###############################################################################
"""
N = 2^14
#sigs = sum(sigs, dims=2)
sigt = [3,3,3]
sigs = [1.5,1.5,1.5]

qmc_data = qmc_init(N, Nx, na2, s, sigs, sigt)
phi_avg, phi_edge, dphi, J_avg, J_edge = qmc_sweep(qmc_data)
midpoints = qmc_data.midpoints
edges = qmc_data.edges

G = size(sigt)[1] # number of groups
Q = 1.0 # source strength
flux = inv(Diagonal(sigt[:,1] - sigs[:,1]))*Q # returns diagonal matrix
flux = sum(sum(flux,dims=2))*ones(Nx) # sum across groups

figure(2)
title("Average Scalar Flux")
plot(midpoints,phi_avg, label="QMC")
plot(midpoints,flux,label="Infinite Medium Sol")
ylabel("cell averaged flux")
xlabel("midpoints")
legend()
"""

###############################################################################
#### Single Function Loop
###############################################################################
"""
N = [2^8, 2^11, 2^14, 2^17]
figure(1, figsize = (20,4))
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

end
"""
