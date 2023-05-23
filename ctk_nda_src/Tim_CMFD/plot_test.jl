function plot_test(N, Nx, Nc)
tau=5.0; dx=tau/Nc; x=.5*dx:dx:tau-.5*dx; s=1.0;
z=dx:dx:tau-dx;
gm_out=cmfd_krylov(N, Nx, Nc)
#dsola=zeros(Nc)
#
# cell centered flux
#
sol=gm_out.sol
J_c = gm_out.J_c
figure(1)
(dphidx, Dhat) = cmfd_parts(sol, J_c, dx)
#
# You can see why we need large N in the plots
#
plot(x,sol,"k-", x, dphidx, "k--", x, J_c, "k.", x, Dhat, "k-.")
title("Flux data: N = $N, Nx = $Nx, Nc = $Nc")
legend([L" \phi", L"d \phi/dx", L"J_c", L"\hat{D}"])
figure(2)
(phi_ct, J_ct) = sn_ok(Nx, Nc, s)
plot(x, J_c,"k-", x, J_ct, "k--")
println("Flux error = ",norm(phi_ct - sol,Inf),"  ", 
        "Current error =", norm(J_ct - J_c, Inf))
phi_sn = zeros(Nc)
s=1.0
phi_sn = sn_verify(Nx, Nc, s)
figure(3)
plot(x,sol-phi_ct)
println(norm(sol-phi_sn)/sqrt(Nc))
end

function cmfd_parts(sol, J_c, dx)
Nc=length(sol)
dsola=zeros(Nc)
#
# derivative of flux at interior cell edges
#
@views dsol=(sol[2:Nc]-sol[1:Nc-1])/dx
#
# Map to cell centers
#
@views dsola[2:Nc-1] = .5*(dsol[1:Nc-2] + dsol[2:Nc-1])
dsola[1]=1.5 * dsol[1]-.5*dsol[2]
dsola[Nc]=1.5 * dsol[Nc-1]-.5*dsol[Nc-2]
#
# Now do the currents and D-hat
#
Dhat = ((dsola/3.0) + J_c)./sol
return(dsola, Dhat)
end

function sn_ok(Nx, Nc, s)
sn_data=sn_init(4096*4, 22, s)
kout=krylov_iteration(sn_data, s, 1.e-8)
phi_e = ctr2edge(kout.sol)
(phi_f, J_f) = sn_fj(kout.sol, sn_data)
println(norm(phi_f-kout.sol,Inf))
phi_c = zeros(Nc)
J_c = zeros(Nc)
ftoc(phi_f, phi_c)
ftoc(J_f, J_c)
return (phi_c, J_c)
end

function sn_verify(Nx, Nc, s)
sn_data=sn_init(4096*4, 22, s)
kout=krylov_iteration(sn_data, s, 1.e-8)
phi_c = zeros(Nc)
ftoc(kout.sol,phi_c)
return phi_c
end
