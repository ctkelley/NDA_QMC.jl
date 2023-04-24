function avplot(N=4096, Nx=512, Nc=64, na2=11; s=1.0, tol=1.e-10)
kout=qmc_krylov(N, Nx, na2; s=s, tol=tol)
uf=kout.sol
uc=zeros(Nc)
ftoc(uf,uc)
plot(uc)
return uc
end
