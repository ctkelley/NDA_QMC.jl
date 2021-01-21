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
#
# precomputed data
#
(angles, weights) = sn_angles(na2)
sn_data=sn_init(nx,na2,s,angles,weights,method)
#
# Problem data
#
c=sn_data.c
psi_right=sn_data.psi_right
psi_left=sn_data.psi_left
psi=sn_data.psi
#
# Source iteration
#
itt=0;
delflux=1;
phi=zeros(nx,)
flux=zeros(nx,)
ithist=[]
while itt < 200 && delflux > 1.e-8
   flux=flux_map!(flux,sn_data)
   delflux=norm(flux-phi,Inf); itt=itt+1;
   push!(ithist,delflux)
   phi.=flux;
end
#
# Tabulate the exit distributions to check results.
#
tout=sn_tabulate(s,nx,flux)
return(left=tout.left, right=tout.right, history=ithist)
end
