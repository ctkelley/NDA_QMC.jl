function result_quality(Nc = 64)
na2=11
Nx=512
dx=5.0/Float64(Nc)
x=.5*dx:dx:5.0-.5*dx
sol=zeros(Nc,3)
N=[1024, 2048, 4096]
#N *= 4
for in=1:3
nout=cmfd_krylov(N[in], Nx, Nc, na2)
sol[:,in]=nout.sol
end
plot(x,sol[:,1],"k-",x,sol[:,2],"k--",x,sol[:,3],"k-.")
N1=N[1]; N2=N[2]; N3=N[3];
legend(["N=$N1", "N=$N2", "N=$N3"])
title("Converged fluxes: Nc=$Nc; Nx=$Nx; s=1.0")
xlabel("x")
ylabel(L"\phi")
end

