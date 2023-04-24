function plot_test(N, Nx, Nc)
tau=5.0; dx=tau/Nc; x=.5*dx:dx:tau-.5*dx;
z=dx:dx:tau-dx;
gm_out=cmfd_krylov(N, Nx, Nc)
dsola=zeros(Nc)
#
# cell centered flux
#
sol=gm_out.sol
#
# derivative of flux at interior cell edges
#
@views dsol=(sol[2:Nc]-sol[1:Nc-1])/dx
#
# Map to cell centers
#
@views dsola[2:Nc-1] = .5*(dsol[1:Nc-2] + dsol[2:Nc-1])
dsola[1]=2.0 * dsol[1]-dsol[2]
dsola[Nc]=2.0 * dsol[Nc-1]-dsol[Nc-2]
#
# Now do the currents and D-hat
#
J_c = gm_out.J_c
Dhat = ((dsola/3.0) + J_c)./sol
#
# You can see why we need large N in the plots
#
plot(x,sol,"k-", x, dsola, "k--", x, J_c, "k.", x, Dhat, "k-.")
title("Flux data: N = $N, Nx = $Nx, Nc = $Nc")
legend([L" \phi", L"d \phi/dx", L"J_c", L"\hat{D}"])
end
