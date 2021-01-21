"""
 source_iteration(na2=20, s=Inf)
 Source iteration example script for transport equation.

 This is one of the test cases from

 Radiative transfer in finite inhomogeneous plane-parallel atmospheres
 by Garcia and Siewert
 JQSRT (27), 1982 pp 141-148.

"""
function source_iteration(na2=20, s=Inf)
nx=2001
method=:new
(angles, weights) = sn_angles(na2)
sn_data=sn_init(nx,na2,s,angles,method)
#tau=5.0; dx=tau/(nx-1);
dx=sn_data.dx
na=floor(Int,na2/2)
#
# Problem data
#
#x=collect(0:dx:tau);
#c=exp.(-x/s);
c=sn_data.c
psi_right=sn_data.psi_right
psi_left=sn_data.psi_left
#psi_left=ones(na,);
#psi_right=zeros(na,);
source=zeros(nx,);
#ptmp=zeros(na,)
ptmp=sn_data.ptmp
sweep_data=(c=c, dx=dx, ptmp=ptmp)
psi=sn_data.psi
#psi=zeros(nx,na2)
#
# Source iteration
#
itt=0;
delflux=1;
phi=zeros(nx,)
flux=zeros(nx,)
ithist=[]
while itt < 200 && delflux > 1.e-8
   psi = transport_sweep!(psi, phi, source, sn_data);
   if method == :old
   flux.=psi*weights;
   else
   flux .= (weights'*psi)'
   end
   delflux=norm(flux-phi,Inf); itt=itt+1;
   push!(ithist,delflux)
   phi.=flux;
end
#
# Tabulate the exit distributions to check results.
#
if true
tout=sn_tabulate(s,nx,flux)
return(left=tout.left, right=tout.right, history=ithist)
else
angleout=[-.05; collect(-.1:-.1:-1.0); .05; collect(.1:.1:1.0)]
na2=length(angleout)
na=floor(Int,na2/2)
sn_data=sn_init(nx, na2, s, angleout, method)
psi=sn_data.psi
psi = transport_sweep!(psi, phi, source, sn_data);
if method == :old
return (left=psi[1,1:na], right=psi[nx,na+1:na2], history=ithist)
else
return (left=psi[1:na,1], right=psi[na+1:na2,nx], history=ithist)
end
end
end
