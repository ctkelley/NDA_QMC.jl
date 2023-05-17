function phi_c = sn_verify(phi_in, cmfd_data)
qmc_data=cmfd_data.qmc_data
phi_f = cmfd_data.phi_f
Nx=length(phi_f)
J_c = cmfd_data.J_c
Nc=length(phi_c)
ctof(phi_c, phi_f)
s=qmc_data.c
sn_data=sn_init(Nx, 22, s; phiedge=false)
kout=krylov_iteration(sn_data,s,tol=1.e-8; onlygmres=true)
phi_out=kout.sol
return phi_out
end
