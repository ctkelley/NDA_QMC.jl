using LinearAlgebra
include("qmc_sweep_2.jl")
"""
 source_iteration(sn_data,s,tol=1.e-8)
 Source iteration example script for transport equation.

 This is one of the test cases from

 Radiative transfer in finite inhomogeneous plane-parallel atmospheres
 by Garcia and Siewert
 JQSRT (27), 1982 pp 141-148.

"""


function qmc_source_iteration(s, qmc_data, tol=1.e-8)
    #
    # precomputed data
    #
    N = qmc_data.N
    Nx = qmc_data.Nx
    itt = 0
    delflux = 1
    phi_avg = zeros(Nx)
    phi_avg_old = zeros(Nx)
    reshist = []

    while itt < 200 && delflux > tol
        phi_avg, phi_edge, J_avg, J_edge, exit_right_bins, exit_left_bins = qmc_sweep(phi_avg,qmc_data)
        delflux = norm(phi_avg - phi_avg_old, Inf)
        itt += 1
        push!(reshist, delflux)
        phi_avg_old .= phi_avg
        #println("**********************")
        #println("Iteration: ", itt," change = ",delflux)
    end

    return (phi_avg = phi_avg,
            phi_edge = phi_edge,
            J_avg = J_avg,
            J_edge = J_edge,
            psi_right = exit_right_bins,
            psi_left = exit_left_bins,
            history = reshist,
            itt = itt)
end
