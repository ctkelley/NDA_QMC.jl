""" 
Fcmfd!(FS, phi, cmfd_data) 
Nonlinear equation form of NDA equations 
""" 
function Fcmfd!(FS, phi,cmfd_data) 
FS.=phi-cmfd_fixed(phi,cmfd_data) 
return FS 
end

function cmfd_nsoli(Nc=64; s=1.0)
FS=zeros(Nc,)
FPS=zeros(Nc,20)
phi0=zeros(Nc,)
N=4096*4; Nx=1024; na2=11;
cmfd_data=cmfd_init(N, Nx, Nc, na2, s);
#
# Fix the initial iterate if you want decent results.
#
phi0.=cmfd_fixed(phi0,cmfd_data)
nout=nsoli(Fcmfd!, phi0, FS, FPS; eta=.1, fixedeta=false,
              rtol=1.e-10, pdata=cmfd_data)
#noutb=nsoli(Fcmfd!, phi0, FS, FPS; eta=.1, fixedeta=false,
#              rtol=1.e-10, lsolver="bicgstab", pdata=cmfd_data)
#return (nout=nout, noutb=noutb)
return nout
end


"""
cmfd_fixed(phi_in, cmfd_data)
Nonlinear fixed point map for NDA

Hardwired Sigma_t = 1; Sigma_s = c

Input:  phi, the low-order flux
        sn_data, data for transport equation and grids

Output: nonlinear nda fixed point map applied to phi
"""
function cmfd_fixed(phi_in, cmfd_data)
#
# Solve high order problem
#
qmc_data=cmfd_data.qmc_data
nc=length(phi_in); dx=5.0/nc
#
# we get the cell average flux, so convert to cell edges
#
phi=copy(phi_in)
phic=copy(phi_in)
(phi, J) = cmfd_sweep(phic, cmfd_data)
#(phi, J) = sn_sweep(phic, cmfd_data)
#bcl=phi[1]; bcr=phi[nc];
bcl=1.5*phi[1]-.5*phi[2]
bcr=1.5*phi[nc]-.5*phi[nc-1]
#
# Build low order problem
#
(dphidx, dhat) = cmfd_parts(phi, J, dx)
#
# Map dhat to cell edges
#
dhat2=ctr2edge(dhat)
Dhat = Diagonal(dhat2)
#
L2=cmfd_data.nda_data.L2
D1=cmfd_data.nda_data.D1
LT = L2 + (D1*Dhat)
# Fix LT to respect the boundary conditions
LT[1,1]=1.0; LT[1,2:nc+1].=0.0;
LT[nc+1,nc+1]=1.0; LT[nc+1,1:nc].=0.0;
#PT=D1*(Dhat*AV)
rhs=zeros(nc+1); rhs[1]=bcl; rhs[nc+1]=bcr;
#
# Solve low order problem for cell-edge flux
#
phi_e=LT\rhs;
#
# Map to cell centers
#
phiout=copy(phi)
@views phiout .= (phi_e[1:nc] + phi_e[2:nc+1])*.5
return phiout
end

function nda_cmfd_init(nx, nc, na2, s)
sn_data=sn_init(nc+1, na2, s)
dx=sn_data.dx
c=sn_data.c
nout=cmfdD2(nc+1,dx,c)
L2=nout.L2
D1=nout.D1
return (sn_data=sn_data, L2=L2, D1=D1)
end

"""
cmfdD2(n,dx,c)

returns -d^2/dx^2 and d/dx on [0,1] zero BC
"""
function cmfdD2(n,dx,c)
    d = 2.0 * ones(n)
    sup = -ones(n - 1)
    slo = -ones(n - 1)
    D2 = Tridiagonal(slo, d, sup)
#
# Fixing Sigma_t = 1 here
#
    D2 = D2 / (3.0 *dx * dx)
#
# Now to capture the boundary conditions
#
#    D2[1,1]=1.0; D2[1,2]=0.0
#    D2[n,n]=1.0; D2[n,n-1]=0.0
#
# First derivative operator
# So it's u' = (U_{i+1} - U_{i-1})/(2h) for 2 <= i <= N-1
# and 0 for i=1 and i=N
# and D1*U = u'
#
    D1s=spdiagm(n, n, -1 => -ones(n-1,), 1 => ones(n-1,))
    D1s ./= (2.0*dx)
#    D1s=spdiagm(n, n, 0 => -ones(n,), 1 => ones(n-1,))
#    D1s ./= dx
#
#   Map to put values at interior points into full vector padded with zeros.
#
    D0=Diagonal(1.0 .- c)
#    D0[1,1]=0.0
#    D0[n,n]=0.0
    L2=D2 + D0
    return (L2 = L2, D1 = D1s)
end

function ctr2edge(cval)
N=length(cval)
eval=zeros(N+1)
eval[2:N]=.5*(cval[1:N-1]+cval[2:N])
eval[1]=2.0*cval[1]-eval[2]
eval[N+1]=2.0*cval[N]-eval[N-1]
#tst1=1.5*cval[1]-.5*cval[2]
#tstN=1.5*cval[N]-.5*cval[N-1]
#println(tst1-eval[1],"  ",tstN-eval[N+1])
return eval
end

function sn_sweep(phi_in, cmfd_data)
sn_data=cmfd_data.sn_data
Nx = sn_data.nx
phi_cf=zeros(Nx-1)
J_cf = copy(phi_cf)
phi_ef=zeros(Nx)
phi_out = copy(phi_in)
J_out = copy(phi_in)
phi_cf .= ctof(phi_in, phi_cf)
phi_ef .= ctr2edge(phi_cf)
(phi_f, J_f) = sn_fj(phi_ef, sn_data)
phi_cf .= (phi_f[1:Nx-1] + phi_f[2:Nx])*.5
J_cf .= (J_f[1:Nx-1] + J_f[2:Nx])*.5
phi_out = ftoc(phi_cf,phi_out)
J_out = ftoc(J_cf, J_out)
return (phi_sn = phi_out, J_sn = J_out)
end

function sn_fj(phi, sn_data)
    J = copy(phi)
    flux=copy(phi)
    psi = sn_data.psi
    psi = transport_sweep!(psi, phi, sn_data)
    #
    # Take the 0th moment to get the flux.
    #
    weights = sn_data.weights
    angles = sn_data.angles
    flux .= (weights' * psi)'
    mom1 = weights.*angles
    J .= (mom1' * psi)'
    return (phi_f = flux, J_f = J)
end

