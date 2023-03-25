function qmc_test(N=10^4, Nx=100, tol=1.e-3, na2=11; s=1.0)
#
qmc_data=qmc_init(N, Nx, na2, s)
nda_data=nda_init(Nx+1, 2*na2, s)
qmc_nda_data=(qmc_data=qmc_data, nda_data=nda_data)
sn_data=nda_data.sn_data
#
# Initial iterate from two deterministic transport sweeps
#
phiD=zeros(Nx+1,)
phiD_ave=zeros(Nx,)
PhiD=flux_map!(phiD,sn_data)
PhiD=flux_map!(phiD,sn_data)
phiD_ave=.5*(phiD[1:Nx]+phiD[2:Nx+1])
#
# Now compare all the stuff. Start with QMC
#
phin=copy(phiD_ave)
qout=qmc_sweep(phin,qmc_data)
phiQ=qout.phi_avg
phiQe=zeros(Nx+1,)
dphiQ=zeros(Nx,1)
JQ=qout.J_avg
avg2edge!(phiQe, dphiQ, phiQ, qmc_nda_data)
#
# And now the deterministic stuff 
#
Dout=all_out(phiD,sn_data)
JD=Dout.JD
println(norm(phiD_ave))
phiD_ave.=Dout.phiD_ave
println(norm(phiD_ave))
dphiD=Dout.dphiD
#
# What's different
#
delphi=norm(phiD_ave-phiQ,Inf)/norm(phiD_ave,Inf)
figure(1)
plot(phiD_ave-phiQ)
delcur=norm(JD-JQ,Inf)/norm(JD,Inf)
figure(2)
plot(JD-JQ)
delder=norm(dphiD-dphiQ,Inf)
println(delphi,"  ",delcur)
end

function all_out(phiD, sn_data)
# Returns cell edge stuff
psi=sn_data.psi
println(norm(psi,Inf))
psi=transport_sweep!(psi, phiD, sn_data)
println(norm(psi,Inf))
weights=sn_data.weights;
angles=sn_data.angles;
nx=sn_data.nx;
dx=sn_data.dx;
hoflux=zeros(nx,);
hocurrent=zeros(nx,);
pa=weights.*angles;
hoflux=(weights' * psi)';
hocurrent=(pa' * psi)';
dflux=(hoflux[2:nx]-hoflux[1:nx-1])
dflux ./= dx
phiav=.5*(hoflux[2:nx] .+ hoflux[1:nx-1])
jav=.5*(hocurrent[2:nx] .+ hocurrent[1:nx-1])
return (JD=jav, phiD_ave=phiav, dphiD=dflux)
end

function avg2edge!(phi_edge, dphi, phi_avg,qmc_nda_data)
sn_data=qmc_nda_data.nda_data.sn_data
dx=sn_data.dx
Nx=length(phi_avg)
m=length(phi_edge)
Nx==length(dphi) || error("dphi lives on cell centers")
m == Nx+1 || error("edges and averages fail to match")
phi_edge[1]=3.0*phi_avg[1] - 2.0*phi_avg[2]
for ie=2:Nx
phi_edge[ie] = .5*(phi_avg[ie-1]+phi_avg[ie])
end
phi_edge[Nx+1]=3.0*phi_avg[Nx] - 2.0*phi_avg[Nx-1]
for id=1:Nx
dphi[id]=(phi_edge[id+1]-phi_edge[id])/dx
end
return (phi_edge=phi_edge, dphi=dphi)
end

function nda_qmc_fixed(phi,qmc_nda_data)
qmc_data=qmc_nda_data.qmc_data
nda_data=qmc_nda_data.nda_data
Nx=length(phi)
phi_edge=zeros(Nx+1,)
dphi=zeros(Nx,)
qout=qmc_sweep(phi,qmc_data)
testmode=true
if testmode
phi_edgep=zeros(Nx+1,)
dphip=zeros(Nx,)
avg2edge!(phi_edgep, dphip, qout.phi_avg, qmc_nda_data)
sn_data=nda_data.sn_data
psi=sn_data.psi
psi=transport_sweep!(psi,phi_edgep,sn_data)
dhatp=zeros(Nx+1,)
(dhatp, bcl1, bcr1)=dhateval!(dhatp, psi, sn_data)
end
avg2edge!(phi_edge, dphi, qout.phi_avg,qmc_nda_data)
dhat=(qout.J_avg + (dphi./3.0))./qout.phi_avg
#println(norm(dhat-dhatp,Inf),"  ",norm(dphi-dphip))
#
# Something's busted here. My normal code has phi at edges. The
# QMC code does centers. I tried Sam's QMC edge values and the
# Taylor's theorem way and nothing seems right.
#
Dhat=Diagonal(dhat)
Nx=length(qout.phi_avg)
L2=nda_data.L2
D1=nda_data.D1
AV=nda_data.AV
LT = L2 + (D1*(Dhat * AV))
# Fix LT to respect the boundary conditions
LT[1,1]=1.0; LT[1,2:nx].=0.0;
LT[nx,nx]=1.0; LT[nx,1:nx-1].=0.0;
bcl=phi_edge[1];
bcr=phi_edge[Nx+1];
#bcl=qout.phi_avg[1];
#bcr=qout.phi_avg[Nx];
rhs=zeros(Nx+1,)
rhs[1]=bcl
rhs[Nx+1]=bcr
if testmode
println(bcl-bcl1,"  ",bcr-bcr1)
end
phiout=LT\rhs
phiave=.5*(phiout[1:Nx]+phiout[2:Nx+1])
#phiout=qout.phi_avg
end


