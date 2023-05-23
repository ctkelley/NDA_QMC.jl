function cmfd_si(N=4096*4, Nx=512, Nc=64, na2=11; s=1.0, tol=1.e-10)
outputok=true
phi_c = zeros(Nc)
phi_c0 = zeros(Nc)
cmfd_data=cmfd_init(N, Nx, Nc, na2, s)
linop_data=linop_init(Nc, cmfd_data)
reshist=Float64[]
resid=1.0
for iz=1:20
   (phi_c, J_c) = cmfd_sweep(phi_c, cmfd_data)
   resid=norm(phi_c - phi_c0, Inf)
   push!(reshist, resid)
   phi_c0 .= phi_c
end
if outputok == true
lin_resid=cmfd_matvec(phi_c, linop_data) - linop_data.frhs
println(norm(lin_resid,Inf))
plot(phi_c)
end
return (sol=phi_c, reshist=reshist)
end


"""
cmfd_krylov(N=4096*4, Nx=512, Nc=64, na2=11; s=1.0, tol=1.e-10)
Solves the linear system you get from CMFD with GMRES.
"""
function cmfd_krylov(N=4096*4, Nx=512, Nc=64, na2=11; s=1.0, tol=1.e-10)
phi_c = zeros(Nc)
phi_c0 = zeros(Nc)
cmfd_data=cmfd_init(N, Nx, Nc, na2, s)
linop_data=linop_init(Nc, cmfd_data)
b=linop_data.frhs
V=zeros(Nc,20)
eta=tol
#
# kl_gmres and kl_bicgstab solve the problem.
#
gout = kl_gmres(phi_c, b, cmfd_matvec, V, tol; pdata = linop_data)
sol=gout.sol
reshistg=gout.reshist
(phi, J_c) = cmfd_sweep(sol, cmfd_data)
return(sol=sol, J_c=J_c, reshistg=reshistg)
end

"""
cmfd_sweep(phi_c, cmfd_data)

Returns coarse mesh cell average flux and current
"""
function cmfd_sweep(phi_c, cmfd_data)
#
# Map coarse mesh flux to fine mesh
#
phi_f = cmfd_data.phi_f
ctof(phi_c,phi_f)
#
# Do Sam's QMC sweep
#
qmc_data=cmfd_data.qmc_data
#
# Is this the thing to do?
# qmc_out.J_avg .*= 0.0
#
qmc_out=qmc_sweep(phi_f, qmc_data)
phi_f .= qmc_out.phi_avg
#
# Map output flux and current to coarse mesh
#
phi_c = ftoc(phi_f, phi_c)
J_f = qmc_out.J_avg
J_c = cmfd_data.J_c
J_c = ftoc(J_f, J_c)
return (phi_c, J_c)
end

"""
cmfd_matvec(phi,linop_data)
Matvec for CMFD-QMC
Does a QMC transport sweep with flux phi, subtracts off the zero bc solution
to obtain the linear integral operator. Then subtract that from the identity.
"""
function cmfd_matvec(phic,linop_data)
phi_copy=linop_data.phi_copy
phi_copy .= phic
frhs=linop_data.frhs
cmfd_data=linop_data.cmfd_data
(phi_c, J_c) =cmfd_sweep(phi_copy,cmfd_data)
Mprod=phic-(phi_c - frhs)
return Mprod
end


function cmfd_init(N, Nx, Nc, na2, s)
phi_f = zeros(Nx)
J_c = zeros(Nc)
sn_data=sn_init(Nx+1, 2*na2, s)
qmc_data = qmc_init(N, Nx, na2, s);
nda_data=nda_cmfd_init(Nx, Nc, 2*na2, s)
return(phi_f=phi_f, J_c = J_c, qmc_data=qmc_data, sn_data=nda_data.sn_data, 
       nda_data=nda_data)
end

function linop_init(Nc, cmfd_data)
nullin=zeros(Nc)
phi_copy=zeros(Nc)
(frhs, Jrhs)=cmfd_sweep(nullin, cmfd_data)
return (frhs = frhs, cmfd_data = cmfd_data, phi_copy = phi_copy)
end
