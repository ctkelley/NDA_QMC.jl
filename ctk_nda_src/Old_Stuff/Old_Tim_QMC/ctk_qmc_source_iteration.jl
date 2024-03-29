#using LinearAlgebra
#include("qmc_sweep.jl")
"""
 qmc_source_iteration(s, qmc_data,tol=1.e-4)
 quasi-monte carlo source iteration example script for transport equation.

 This is one of the test cases from

 Radiative transfer in finite inhomogeneous plane-parallel atmospheres
 by Garcia and Siewert
 JQSRT (27), 1982 pp 141-148.

 Inputs:
 -------
     N:
         Number of internal source particles
     Nx:
         Number of scalar flux and current tallying cells
     na2:
         Number of exit angular flux tally bins
     s:
         Scattering cross section varies with exp(-x/s)
 Outputs:
 --------
     phi_avg:
         Array of length Nx, average scalar flux across each cell
     phi_edge:
         Array of length Nx+1, scalar flux at cell boundaries
     J_avg:
         Array of length Nx, average cell current
     J_edge:
         Array of length Nx+1, current at cell boundaries
     psi_right:
         Array of length nbins, angular flux distribution at right boundary
     psi_left:
         Array of length nbins, angular flux distribution at left boundary
     history:
         Infinity norm of change in phi_avg with each iteration
     itt:
         Number of iterations
"""


function qmc_source_iteration(s, qmc_data, tol=1.e-4)
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
    #globalize variables
    phi_edge = J_avg = J_edge = exit_right_bins = exit_left_bins = 0

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
