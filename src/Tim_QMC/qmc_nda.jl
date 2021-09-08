function qmc_nda(N=10^4, Nx=100, tol=1.e-3, na2=11; s=1.0, 
       nda=true, tabulate=false)
#
qmc_data=qmc_init(N, Nx, na2, s)
nda_data=nda_init(Nx+1, 2*na2, s)
qmc_nda_data=(qmc_data=qmc_data, nda_data=nda_data)
#
# Initial iterate from one transport sweep
#
phin=zeros(Nx,)
phiout=ones(Nx,)
itmax=30
itc=0
ithist=[]
phidiff=norm(phin-phiout,Inf)
testmode=~nda
while (phidiff > tol) && itc < itmax
if testmode
qout=qmc_sweep(phin,qmc_data)
phiout=qout.phi_avg
else
phiout=nda_qmc_fixed(phin,qmc_nda_data)
end
phidiff=norm(phin-phiout,Inf)
phin .= phiout
push!(ithist,phidiff)
itc+=1
end
if tabulate
if nda
println("QMC-NDA Exit Distributions")
else
println("QMC Exit Distributions")
end
sn_tabulate(s, Nx, phiout; phiedge=false)
end
return(phiout=phiout, ithist=ithist)
end

function avg2edge!(phi_edge, dphi, phi_avg,qmc_nda_data)
sn_data=qmc_nda_data.nda_data.sn_data
dx=sn_data.dx
Nx=length(phi_avg)
m=length(phi_edge)
Nx==length(dphi) || error("dphi lives on cell centers")
m == Nx+1 || error("edges and averages fail to match")
phi_edge[1]=1.5*phi_avg[1] - .5*phi_avg[2]
for ie=2:Nx
phi_edge[ie] = .5*(phi_avg[ie-1]+phi_avg[ie])
end
phi_edge[Nx+1]=1.5*phi_avg[Nx] - .5*phi_avg[Nx-1]
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
bcl=phi_edge[1];
bcr=phi_edge[Nx+1];
#bcl=qout.phi_avg[1];
#bcr=qout.phi_avg[Nx];
rhs=zeros(Nx+1,)
rhs[1]=bcl
rhs[Nx+1]=bcr
phiout=LT\rhs
phiave=.5*(phiout[1:Nx]+phiout[2:Nx+1])
#phiout=qout.phi_avg
end

