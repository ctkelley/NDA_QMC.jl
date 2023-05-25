function test_cmfd(N=4096*4, Nx=1024, Nc=64; s=1.0, tau=5.0)
#N=4096*4; Nx=1024; 
na2=11; 
#N=4096*2; Nx=1024; na2=11;
dx=tau/Nc; x=.5*dx:dx:tau-.5*dx
cmfd_data=cmfd_init(N, Nx, Nc, na2, s; tau=tau);
phi0=zeros(Nc);
phic=copy(phi0)
reshist=Float64[]
plotit=false
itc=0
resid=1.0
while (itc < 200) && (resid > 1.e-11)
itc+=1
    phip=cmfd_fixed(phic, cmfd_data)
resid=norm(phip-phic)
push!(reshist,resid)
if plotit && itc <=10
println(norm(phip-phic,Inf),"   ",phip[end],"  ",phip[end-1],"  ",phip[end-2])
plot(x,phip)
end
    phic .= phip
end
return (solution=phic, history=reshist)
end
