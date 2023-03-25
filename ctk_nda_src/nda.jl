"""
Fnda!(FS, phi, nda_data)

Nonlinear equation form of NDA equations
"""
function Fnda!(FS, phi,nda_data)
FS.=phi-nda_fixed(phi,nda_data)
return FS
end


"""
nda_fixed(phi, nda_data)
Nonlinear fixed point map for NDA

Hardwired Sigma_t = 1; Sigma_s = c

Input:  phi, the low-order flux
        sn_data, data for transport equation and grids

Output: nonlinear nda fixed point map applied to phi
"""
function nda_fixed(phi, nda_data)
#
# Solve high order problem
#
sn_data=nda_data.sn_data;
weights=sn_data.weights;
angles=sn_data.angles;
nx=sn_data.nx;
dx=sn_data.dx;
na=length(angles);
psi=sn_data.psi;
psi=transport_sweep!(psi, phi, sn_data);
dhat=zeros(nx-1,)
(dhat, bcl, bcr)  = dhateval!(dhat, psi, sn_data, phi);
Dhat = Diagonal(dhat)
L2=nda_data.L2
D1=nda_data.D1
AV=nda_data.AV
LT = L2 + (D1*(Dhat * AV))
# Fix LT to respect the boundary conditions
LT[1,1]=1.0; LT[1,2:nx].=0.0;
LT[nx,nx]=1.0; LT[nx,1:nx-1].=0.0;
#PT=D1*(Dhat*AV)
rhs=zeros(nx,); rhs[1]=bcl; rhs[nx]=bcr;
phiout=LT\rhs;
return phiout
end

function dhateval!(dhat, psi, sn_data, phi)
weights=sn_data.weights;
angles=sn_data.angles;
nx=sn_data.nx;
dx=sn_data.dx;
hoflux=zeros(nx,);
hocurrent=zeros(nx,);
pa=weights.*angles;
hoflux=(weights' * psi)';
bcl=hoflux[1]; bcr=hoflux[nx];
#bcl=(hoflux[1]+hoflux[2])*.5
#bcr=(hoflux[nx]+hoflux[nx-1])*.5
hocurrent=(pa' * psi)';
#println("flux error = ",norm(hoflux-phi,Inf))
dflux=zeros(nx-1,)
dhat=zeros(nx-1,)
phiav=zeros(nx-1,)
jav=zeros(nx-1,)
dflux.=(hoflux[2:nx]-hoflux[1:nx-1])
dflux./= (3.0*dx)
phiav.=.5*(hoflux[2:nx] .+ hoflux[1:nx-1])
jav.=.5*(hocurrent[2:nx] .+ hocurrent[1:nx-1])
dhat.=(jav + dflux)./phiav
return (dhat, bcl, bcr)
end


function nda_init(nx, na2, s)
sn_data=sn_init(nx, na2, s)
dx=sn_data.dx
c=sn_data.c
nout=ndaD2(nx,dx,c)
L2=nout.L2
D1=nout.D1
AV=nout.AV
return (sn_data=sn_data, L2=L2, D1=D1, AV=AV)
end

"""
ndaD2(n,dx,c)

returns -d^2/dx^2 and d/dx on [0,1] zero BC
"""
function ndaD2(n,dx,c)
    d = 2.0 * ones(n)
    sup = -ones(n - 1)
    slo = -ones(n - 1)
    D2 = Tridiagonal(slo, d, sup)
#
# Fixing Sigma_t = 1 here.j
#
    D2 = D2 / (3.0 *dx * dx)
#
# Now to capture the boundary conditions
#
    D2[1,1]=1.0; D2[1,2]=0.0
    D2[n,n]=1.0; D2[n,n-1]=0.0
#
# First derivative operator mapping cell averages to interior cell edges
# So it's u' = (U_{i+1} - U_i)/h if $U_i = (u_{i+1} + u_i)/2
# and D1*U = u'
# The size of D1 is (n-2) x (n-1)
#
    D1s=spdiagm(n-2, n-1, 0 => -ones(n-2,), 1 => ones(n-2,))
    D1s ./= dx
#
# Map from cell edge values to cell average values
#
    AV = spdiagm(n-1, n, 0 => ones(n-1,), 1 => ones(n-1,))
    AV .*= .5
#
#   Map to put values at interior points into full vector padded with zeros.
#
    PadV = spdiagm(n, n-2, -1 => ones(n-2,))
    D1=PadV*D1s
    D0=Diagonal(1.0 .- c)
    L2=D2 + D0
    return (L2 = L2, D2=D2, D1 = D1, AV=AV, PadV=PadV)
end

function testnda(n,k=1)
       sn_data=sn_init(n,10,1.0)
       c=sn_data.c
       dx=1.0/(n-1.0)
       nout=ndaD2(n,dx,c)
       x=collect(0:dx:1.0)
       y=sin.(k*pi*x)
       ya=zeros(n-1,)
       yb=zeros(n-1,)
       AV=nout.AV
       ya .= AV*y
       dy=norm(ya-yb)
       z=k*pi*cos.(k*pi*x)
       D1=nout.D1
#       PadV=nout.PadV
#       p=PadV*(D1*ya)
       p=D1*ya
       derr=zeros(n-2,)
       derr=p[2:n-1]-z[2:n-1]
       D2=nout.D2
       q=D2*y
       d2err=zeros(n,)
       eval=(k*k*pi*pi)/3.0
       d2err[2:n-1]=q[2:n-1]-eval*y[2:n-1]
       L2=nout.L2
       ly=L2*y
       lye=eval*y + (1.0 .- c).*y
       l2err=zeros(n-2,)
       l2err.=ly[2:n-1] - lye[2:n-1]
       println(norm(derr,Inf),"  ",norm(d2err,Inf),"  ",norm(l2err,Inf))
       end


