function ndaerr(Nx=2048, na=40, s=1.0; N=10000)
sn_data=sn_init(Nx, na, s)
kout = krylov_iteration(sn_data,s,1.e-9)
nda_data=nda_init(Nx,na,s)
sn_sol=kout.sol
sn_nda=nda_fixed(sn_sol,nda_data)
ndaerr=sn_sol-sn_nda;
println(norm(ndaerr,Inf))
qmc_data=qmc_init(N, Nx, na, s)
AVE=averop(4,Nx)
qout=qmc_krylov(N, Nx, na)
qsolavg=AVE*(qout.sol)
qsol=zeros(Nx)
qsol[1]=1.5*qsolavg[1] - .5*qsolavg[2]
for ie=2:Nx
qsol[ie] = .5*(qsolavg[ie-1]+qsolavg[ie])
end
qsol[Nx]=1.5*qsolavg[Nx-1] - .5*qsolavg[Nx-2]
qmcsndif=qsol - sn_sol
println(norm(qmcsndif,Inf))
end
