
"""
     qmc_source_iteration(s, qmc_data,tol=1.e-8)
 quasi-monte carlo source iteration example script for transport equation.
"""
function qmc_source_iteration(qmc_data, tol=1.e-5, maxIter=20)
    # precomputed data
    phi_avg = qmc_data.phi_avg
    phi_avg_old = qmc_data.phi_avg
    source = qmc_data.source
    sigs = qmc_data.sigs

    itt = 0
    delflux = 1
    reshist = []
    # initialize global variables
    phi_edge = dphi = J_avg = J_edge = exit_right_bins = exit_left_bins = 0

    while itt < maxIter && delflux > tol
        #q = phi_avg*sigs' + source # multi group data
        q =  phi_avg.*sigs + source # garcia tests and infinite medium problems
        phi_avg, phi_edge, dphi, J_avg, J_edge, exit_right_bins, exit_left_bins = qmc_sweep(q, qmc_data)
        delflux = norm(phi_avg - phi_avg_old, Inf)
        itt += 1
        push!(reshist, delflux)
        phi_avg_old .= phi_avg
        println("**********************")
        println("Iteration: ", itt," change = ",delflux)
    end

    return (phi_avg = phi_avg,
            phi_edge = phi_edge,
            dphi = dphi,
            J_avg = J_avg,
            J_edge = J_edge,
            exit_right_bins = exit_right_bins,
            exit_left_bins = exit_left_bins,
            history = reshist,
            itt = itt)
end
