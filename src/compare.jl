function compare(s = 1.0; qmc=false)
    na2 = 20
    N=1000
    #
    # precomputed data
    #
    if qmc
    N=1000
    nx = 100
    else
    nx = 2001
    sn_data = sn_init(nx, na2, s)
    end
    #
    # do source iteration
    #
    if qmc
    sout = qmc_si(N, nx, na2; s=s, maxit=200)
    fout = sout.sol
    else
    sout = source_iteration(sn_data, s)
    fout = sout.flux
    end
    #
    # and now GMRES
    #
    if qmc
    gout=qmc_gmres(N, nx, na2; s=s)
    else
    gout=gmres_iteration(sn_data,s)
    end
    fgm = gout.sol
    #
    # Make a nice plot
    #
    snhist = sout.history ./ sout.history[1]
    ghist = gout.reshist ./ gout.reshist[1]
    semilogy(ghist, "k-", snhist, "k--")
    legend(["gmres", "source iteration"])
    xlabel("Transport Sweeps")
    ylabel("Relative Residual")
    s == Inf ? strs = L"\infty" : strs = string(s)
    if qmc
    tstring = string("QMC Residual Histories: s= ", strs)
    else
    tstring = string("Residual Histories: s= ", strs)
    end
    title(tstring)
    # sanity check. The difference in solutions should be consistent
    # with the tolerances.
    fgm=gout.sol
    solverdelta=norm(fout - fgm, Inf)
#    (solverdelta < 1.e-6) || error("results not the same")
    (solverdelta < 1.e-6) || println(solverdelta)
    sn_tabulate(s,nx,fgm)
    println("Norm of result difference = ", norm(fout - fgm, Inf))
end
