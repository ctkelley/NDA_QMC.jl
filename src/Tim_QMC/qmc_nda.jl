function qmc_nda(N=10^4, Nx=100, tol=1.e-3, na2=11; s=1.0)
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
while (phidiff > tol) && itc < itmax
qout=qmc_sweep(phin,qmc_data)
phiout=nda_qmc_fixed(phin,qmc_nda_data)
phidiff=norm(phin-phiout,Inf)
phin .= phiout
push!(ithist,phidiff)
itc+=1
end
sn_tabulate(s, Nx, phiout; phiedge=false)
return(phiout=phiout, ithist=ithist)
end


function nda_qmc_fixed(phi,qmc_nda_data)
qmc_data=qmc_nda_data.qmc_data
nda_data=qmc_nda_data.nda_data
Nx=length(phi)
qout=qmc_sweep(phi,qmc_data)
dhat=(qout.J_avg + (qout.dphi./3.0))./qout.phi_avg
#
# Something's busted here. My normal code has phi at edges. The
# QMC code does centers and I need to sort this out.
#
Dhat=Diagonal(dhat)
Nx=length(qout.phi_avg)
L2=nda_data.L2
D1=nda_data.D1
AV=nda_data.AV
LT = L2 + (D1*(Dhat * AV))
bcl=qout.phi_edge[1];
bcr=qout.phi_edge[Nx+1];
#bcl=qout.phi_avg[1];
#bcr=qout.phi_avg[Nx];
rhs=zeros(Nx+1,)
rhs[1]=bcl
rhs[Nx+1]=bcr
phiout=LT\rhs
phiave=.5*(phiout[1:Nx]+phiout[2:Nx+1])
#phiout=qout.phi_avg
end


