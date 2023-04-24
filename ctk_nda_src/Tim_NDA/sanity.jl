function sanity(Nx=2048, na=40, s=1.0, N=10^4)
#
# NDA and NDA-Newton
#
ndaout=nda_iteration(Nx, na, s)
ndakout=nda_nsoli(Nx, na, s)
ndadiff= norm(ndaout.sol - ndakout.nout.solution, Inf)
ndaok = (ndadiff < 1.e-8)
ndaok || println("NDA consistency error = $ndadiff")
#
# SN and SN-Krylov
#
sn_data=sn_init(Nx, na, s)
kout = krylov_iteration(sn_data,s,1.e-9)
sout = source_iteration(sn_data,s,Nx,1.e-9)
sndiff=norm(kout.sol-sout.flux, Inf)
snok= (sndiff < 1.e-8)
snok || println("SN consistency error")
#
# Compare SN and NDA
#
snvsnda=kout.sol-ndaout.sol
snndadiff=norm(snvsnda,Inf)
snndaok=(snndadiff<2.e-6)
snndaok || println("SN vs NDA consistency error")
#
# Compare SN and QMC
#
qout=qmc_krylov(N, Nx, na)
qmdiff=qout.sol-kout.sol
println(norm(qmdiff,2)/sqrt(Float64(Nx)))
#
# Now see if there is any hope for a correct current
#
qmc_data=qmc_init(N, Nx, na, s)
qdout=qmc_sweep(kout.sol,qmc_data)
qediff=qdout.phi_avg-kout.sol
subplot(121)
semilogy(qediff)
println(norm(qediff,2)/sqrt(Float64(Nx)))
(flux, current) = getJ(kout.sol, sn_data)
cdiff=qdout.J_avg - current
println(norm(cdiff,2)/sqrt(Float64(Nx)))
subplot(122)
semilogy(cdiff)
#
return (snndaok && snok && ndaok)
end

function getJ(phi,sn_data)
    psi = sn_data.psi
    psi = transport_sweep!(psi, phi, sn_data)
    #
    # Take the 0th moment to get the flux.
    #
    weights = sn_data.weights
    angles = sn_data.angles
    flux=copy(phi)
    flux .= (weights' * psi)'
    current = copy(phi)
    current .= ((weights.*angles)' * psi)'
    return (flux=flux, current=current)
end
