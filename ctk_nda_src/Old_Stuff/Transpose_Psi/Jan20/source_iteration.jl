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
(angles, weights) = sn_angles(na2)
tau=5.0; dx=tau/(nx-1);
na=floor(Int,na2/2)
#
# Problem data
#
x=collect(0:dx:tau);
c=exp.(-x/s);
psi_right=zeros(na,);
psi_left=ones(na,);
source=zeros(nx,);
ptmp=zeros(na,)
sweep_data=(c=c, dx=dx, ptmp=ptmp)
psi=zeros(nx,na2)
#
# Source iteration
#
itt=0;
delflux=1;
phi=zeros(nx,)
flux=zeros(nx,)
ithist=[]
while itt < 200 && delflux > 1.e-8
   psi = transport_sweep!(psi, phi, psi_left, psi_right, 
              angles, source, sweep_data);
   flux.=psi*weights;
   delflux=norm(flux-phi,Inf); itt=itt+1;
   push!(ithist,delflux)
   phi.=flux;
end
#
# Tabulate the exit distributions to check results.
#
#angleout=[collect(-1.0:.1:-.1); -.05; .05; collect(.1:.1:1.0)]
angleout=[-.05; collect(-.1:-.1:-1.0); .05; collect(.1:.1:1.0)]
na2=length(angleout)
na=floor(Int,na2/2)
psi_left=ones(na,);
psi_right=zeros(na,);
psi = zeros(nx,na2)
ptmp=zeros(na,)
sweep_data=(c=c, dx=dx, ptmp=ptmp)
psi = transport_sweep!(psi, phi, psi_left, psi_right, 
                            angleout, source, sweep_data);
return (left=psi[1,1:na], right=psi[nx,na+1:na2], history=ithist)
end
