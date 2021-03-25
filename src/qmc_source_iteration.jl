"""
 source_iteration(sn_data,s,tol=1.e-8)
 Source iteration example script for transport equation.

 This is one of the test cases from

 Radiative transfer in finite inhomogeneous plane-parallel atmospheres
 by Garcia and Siewert
 JQSRT (27), 1982 pp 141-148.

"""
function source_iteration(N, Nx, na2, s, sn_data, tol=1.e-8)
    #
    # precomputed data
    #
    itt = 0
    delflux = 1
    phi = zeros(Nx)
    flux = zeros(Nx)
    reshist = []
    while itt < 200 && delflux > tol
        flux = flux_map!(flux, sn_data)
        delflux = norm(flux - phi, Inf)
        itt += 1
        push!(reshist, delflux)
        phi .= flux
    end
    #
    # Tabulate the exit distributions to check results.
    #
    return ( flux = flux, history= reshist)
end
