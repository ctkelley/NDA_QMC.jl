function compare(s = 1.0)
    na2 = 20
    nx = 2001
    #
    # precomputed data
    #
    (angles, weights) = sn_angles(na2)
    sn_data = sn_init(nx, na2, s, angles, weights)
    #
    # do source iteration
    sout = source_iteration(sn_data, s)
    fout = sout.flux
if true
    gout=gmres_iteration(sn_data,s)
else
    #
    # Set up the linear system for GMRES
    #
    # Separate the right side from the linear operator
    # First step is apply the fixed point map with the real 
    # boundary conditions and zero flux.
    #
    rhs = zeros(size(fout))
    rhs = flux_map!(rhs, sn_data)
    #
    # Now set the boundary conditions to zero in the
    # precomputed data. That sets up the linear operator.
    #
    sn_data.psi_left .*= 0.0
    sn_data.psi_right .*= 0.0
    #
    # Get organized to call GMRES. Allocate storage for the solution
    # and the Krylov basis
    #
    f0 = zeros(size(fout))
    V = zeros(nx, 20)
    #
    # kl_gmres solves the problem.
    #
    gout = kl_gmres(f0, rhs, sn_matvec, V, 1.e-10; pdata = sn_data)
    #
    #
    # Organize the residual history and make the plot.
    #
end
    fgm = gout.sol
    snhist = sout.history ./ sout.history[1]
    ghist = gout.reshist ./ gout.reshist[1]
    semilogy(ghist, "k-", snhist, "k--")
    legend(["gmres", "source iteration"])
    xlabel("Transport Sweeps")
    ylabel("Relative Residual")
    s == Inf ? strs = L"\infty" : strs = string(s)
    tstring = string("Residual Histories: s= ", strs)
    title(tstring)
    # sanity check. The difference in solutions should be consistent
    # with the tolerances.
    fgm=gout.sol
    fout=
    println(norm(fout - fgm, Inf))
end

"""
sn_matvec(f, sn_data)

Take the source iteration map and extract the
linear operator part. This is the mat-vec for GMRES.
"""
function sn_matvec(f, sn_data)
    kf = copy(f)
    kf = flux_map!(kf, sn_data)
    kf = f - kf
    return kf
end
