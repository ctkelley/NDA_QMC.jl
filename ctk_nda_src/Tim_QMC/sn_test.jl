function sn_solve(Nx=100, tol=1.e-3, na2=22; s=1.0, tabulate=false)
# Nx cells to compare with the QMC code
n=Nx+1
sn_data=sn_init(n, na2, s)
sout=source_iteration_test(sn_data,s)
phiout=sout.flux
end
"""
 source_iteration_test(sn_data,s,tol=1.e-8)
 Source iteration example script for transport equation.

 This is one of the test cases from

 Radiative transfer in finite inhomogeneous plane-parallel atmospheres
 by Garcia and Siewert
 JQSRT (27), 1982 pp 141-148.

"""
function source_iteration_test(sn_data,s,tol=1.e-8)
    nx = sn_data.nx
    #
    # precomputed data
    #
    angles=sn_data.angles
    weihts=sn_data.weights
    itt = 0
    delflux = 1
    phi = zeros(nx)
    flux = zeros(nx)
    reshist = []
    while itt < 200 && delflux > tol
        flux = flux_map!(flux, sn_data)
        delflux = norm(flux - phi, Inf)
        itt = itt + 1
        push!(reshist, delflux)
        phi .= flux
    end
    #
    # Tabulate the exit distributions to check results.
    #
    return ( flux = flux, history= reshist)
end

