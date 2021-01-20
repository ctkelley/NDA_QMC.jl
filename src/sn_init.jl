function sn_init(nx, na2, s, angles, method::Symbol)
tau=5.0; dx=tau/(nx-1);
na=floor(Int,na2/2)
x=collect(0:dx:tau);
c=exp.(-x/s);
psi_right=zeros(na,);
psi_left=ones(na,);
source=zeros(nx,);
ptmp=zeros(na,)
if method == :old
psi=zeros(nx,na2)
else
psi=zeros(na2,nx)
end
return sn_data=(c=c, dx=dx, psi=psi, angles=angles, method=method,
  psi_right=psi_right, psi_left=psi_left, ptmp=ptmp)
end
