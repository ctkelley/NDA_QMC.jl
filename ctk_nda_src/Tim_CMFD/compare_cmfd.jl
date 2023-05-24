function compare_cmfd(N=4096*4, Nx=1024, Nc=64, na2=11; s=1.0, tol=1.e-10)
phi_c = zeros(Nc)
phi_c0 = zeros(Nc)
cmfd_data=cmfd_init(N, Nx, Nc, na2, s)
#
# QMC-GMRES is the baseline 
#
goutgm=qmc_krylov(N, Nx, na2; s=s)
reshistg=goutgm.reshistg
itgmax=length(reshistg)
#
# Newton-GMRES + NDA
#
nout=cmfd_nsoli(Nc; s=s)
fout=nout.stats.ifun+nout.stats.ijac+nout.stats.iarm
#
# count the warm-up qmc-si sweep
#
fout[2]+=1
ftot=cumsum(fout)
reshistn=nout.history
#
# NDA-Picard
#
pout=test_cmfd(Nc; s=s)
reshistp=pout.history
itpic=length(reshistp)
#
# Plot the results
#
semilogy(1:itgmax,reshistg/reshistg[1],"k-",
               1:itpic,reshistp/reshistp[1],"k.-",
               ftot,reshistn/reshistn[1],"k--")
title("N=$N, Nx=$Nx, Nc=$Nc")
xlabel("Transport sweeps")
legend(["QMC-GMRES","NDA-Picard","NDA-Newton"])

end

