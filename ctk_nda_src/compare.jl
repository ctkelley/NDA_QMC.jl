"""
compare(s = 1.0; qmc=false)

Makes the Krylov comparison tables/figures
"""
function compare(s = 1.0; qmc=false)
    na2 = 20
    N=1000
    #
    # precomputed data
    #
    if qmc
    N=1025
    nx = 100
    else
    nx = 2049
    sn_data = sn_init(nx, na2, s)
    end
    #
    # do source iteration
    #
    if qmc
    sout = qmc_si(N, nx, na2; s=s, maxit=200)
    fout = sout.sol
    else
    sout = source_iteration(sn_data, s, nx)
    fout = sout.flux
    end
    #
    # and now GMRES and Bi-CGSTAB
    #
    if qmc
    gout=qmc_krylov(N, nx, na2; s=s)
    else
    gout=krylov_iteration(sn_data,s)
    end
    fgm = gout.sol
    fgmb = gout.solb
    #
    # Make a nice plot
    #
    snhist = sout.history ./ sout.history[1]
    slen=collect(1:length(snhist))
    ghist = gout.reshistg ./ gout.reshistg[1]
    glen=collect(1:length(ghist))
    bhist = gout.reshistb ./ gout.reshistb[1]
    blen=2*collect(1:length(bhist))
    blen[1]=1
    semilogy(glen,ghist, "k-", slen,snhist, "k--", blen, bhist, "k-.")
    legend(["gmres", "source iteration", "bicgstab"])
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
    println("Norm of result differences: GMRES = ", norm(fout - fgm, Inf),
            ",    Bi-CGSTAB = ",norm(fout - fgmb, Inf))
end
