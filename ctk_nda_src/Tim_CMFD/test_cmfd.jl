function test_cmfd(Nc=64)
N=4096*4; Nx=1024; na2=11; s=1.0;
cmfd_data=cmfd_init(N, Nx, Nc, na2, s);
phi0=zeros(Nc);
phic=copy(phi0)
for itc=1:10
    phip=cmfd_fixed(phic, cmfd_data)
println(norm(phip-phic,Inf),"   ",phip[end],"  ",phip[end-1],"  ",phip[end-2])
    phic .= phip
plot(phip)
end
end
