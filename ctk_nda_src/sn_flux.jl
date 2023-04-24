function sn_flux(Nx=4096, Nc=64, na2=22; s=1.0, tol=1.e-10)
sn_data = sn_init(Nx+1, na2, s)
gout=krylov_iteration(sn_data,s)
phi_f=zeros(Nx)
typeof(phi_c)
phi_f .= .5*(gout.sol[2:Nx+1] + gout.sol[1:Nx])
uc=zeros(Nc)
sn_sol_c = ftoc(phi_f, uc)
return uc
end

