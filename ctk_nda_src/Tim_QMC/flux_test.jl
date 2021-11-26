function flux_test()
Nx=100;
na2=11;
s=1.0;
phit=ones(Nx,);
data1=qmc_init(1000,Nx,na2,s);
data2=qmc_init(10000,Nx,na2,s);
phiout1=qmc_sweep(phit, data1);
phiout2=qmc_sweep(phit, data2);
plot(phiout1.phi_edge,"k-");
plot(phiout2.phi_edge,"k--");
legend(["N=1000", "N=10000"])
end
