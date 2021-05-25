function ctk_qmc_nda_test(N,Nx)
s=1.0
na=20
nx=Nx+1
nda_data=nda_init(nx,na,s)
sn_data=nda_data.sn_data
angles=sn_data.angles;
nx=sn_data.nx;
dx=sn_data.dx;
na=length(angles);
psi=sn_data.psi;
phi=zeros(nx,)
phix=copy(phi)
phiout=flux_map!(phix,sn_data)
psix=copy(psi); phiy=copy(phi)
psiout=transport_sweep!(psix, phiy, sn_data; phiedge=true)
dhat=zeros(nx-1,)
classout=fluxdata(psiout,sn_data)
qmc_data=qmc_init(N, Nx, 10, 1.0)
phiq=nda_data.AV*phi
qout=qmc_sweep(phiq, qmc_data)
phip=nda_data.AV*qout.phi_edge
phix=nda_data.AV*phiout
println("L infty absolute errors")
println("Flux edge difference = ", norm(qout.phi_edge-phiout,Inf))
println("Flux average difference = ", norm(phip-phix,Inf))
println("Current edge difference = ",norm(classout.hocurrent-qout.J_edge,Inf))
Jx=nda_data.AV*classout.hocurrent
println("Current average difference = ",norm(Jx-qout.J_avg,Inf))
println("dphi difference  = ", norm(3.0*classout.dflux-qout.dphi,Inf))
# relative l infty errors
println("L infty relative errors")
println("Flux edge difference = ", relerr(qout.phi_edge,phiout,Inf))
println("Flux average difference = ", relerr(phip,phix,Inf))
println("Current edge difference = ",relerr(classout.hocurrent,qout.J_edge,Inf))
println("Current average difference = ",relerr(Jx,qout.J_avg,Inf))
println("dphi difference  = ", relerr(3.0*classout.dflux,qout.dphi,Inf))
# relative l2 errors
println("L2 relative errors")
println("Flux edge difference = ", relerr(qout.phi_edge,phiout))
println("Flux average difference = ", relerr(phip,phix))
println("Current edge difference = ", relerr(classout.hocurrent,qout.J_edge))
println("Current average difference = ",relerr(Jx,qout.J_avg))
println("dphi difference  = ", relerr(3.0*classout.dflux,qout.dphi))
end

function relerr(xf,yf,p=2)
df=xf-yf
relerr=norm(df,p)/norm(yf,p)
return relerr
end

function fluxdata(psi,sn_data)
weights=sn_data.weights;
angles=sn_data.angles;
nx=sn_data.nx;
dx=sn_data.dx;
hoflux=zeros(nx,);
hocurrent=zeros(nx,);
pa=weights.*angles;
hoflux=(weights' * psi)';
bcl=hoflux[1];
bcr=hoflux[nx];
hocurrent=(pa' * psi)';
dflux=zeros(nx-1,)
dhat=zeros(nx-1,)
phiav=zeros(nx-1,)
jav=zeros(nx-1,)
dflux.=(hoflux[2:nx]-hoflux[1:nx-1])
dflux./= (3.0*dx)
phiav.=.5*(hoflux[2:nx] .+ hoflux[1:nx-1])
jav.=.5*(hocurrent[2:nx] .+ hocurrent[1:nx-1])
dhat.=(jav + dflux)./phiav
return (dhat=dhat, bcl=bcl, bcr=bcr, dflux=dflux, 
       hoflux=hoflux, hocurrent=hocurrent)
end


