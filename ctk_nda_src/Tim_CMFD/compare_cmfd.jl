function compare_cmfd(N=4096*4, Nx=1024, Nc=64, na2=11; 
                    tau=5.0, s=1.0, tol=1.e-10, plotme=true)
phi_c = zeros(Nc)
phi_c0 = zeros(Nc)
dx=tau/Nc
x=dx*.5:dx:tau-dx*.5
cmfd_data=cmfd_init(N, Nx, Nc, na2, s; tau=tau)
#
# QMC-GMRES is the baseline 
#
goutgm=qmc_krylov(N, Nx, na2; s=s, tau=tau)
solgmf=goutgm.sol
solgm=zeros(Nc)
solgm.=ftoc(solgmf,solgm)
reshistg=goutgm.reshistg
itgmax=length(reshistg)
#
# Newton-GMRES + NDA
#
nout=cmfd_nsoli(N, Nx, Nc; s=s, tau=tau)
soln=nout.solution
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
pout=test_cmfd(N, Nx, Nc; s=s, tau=tau)
solp=pout.solution
println(norm(soln-solp,Inf),"  ",norm(soln-solgm,Inf))
reshistp=pout.history
itpic=length(reshistp)
#
# Plot the results
#
Itau=Int(tau)
if plotme
figure(1)
plot(x,soln,"k-",x,solgm,"k--",x,solp,"k-.")
figure(2)
semilogy(1:itgmax,reshistg/reshistg[1],"k-",
               1:itpic,reshistp/reshistp[1],"k.-",
               ftot,reshistn/reshistn[1],"k--")
stri=L"\infty"
xxst="N=$N, Nx=$Nx, Nc=$Nc, s="*stri*", Lx = $Itau"
if s==Inf
#title("N=$N, Nx=$Nx, Nc=$Nc, s=, Lx = $Itau")
title(xxst)
else
title("N=$N, Nx=$Nx, Nc=$Nc, s=1.0, Lx = $Itau")
end
xlabel("Transport sweeps")
legend(["QMC-GMRES","NDA-Picard","NDA-Newton"])
end
end

