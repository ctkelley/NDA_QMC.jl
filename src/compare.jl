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
    #
    sout = source_iteration(sn_data, s)
    fout = sout.flux
    #
    # and now GMRES
    #
    gout=gmres_iteration(sn_data,s)
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
    tstring = string("Residual Histories: s= ", strs)
    title(tstring)
    # sanity check. The difference in solutions should be consistent
    # with the tolerances.
    fgm=gout.sol
    println(norm(fout - fgm, Inf))
    #
    # Tabulate the results to compare with the Garcia/Siewert paper
    #
end
