function sanity(Nx=1024, na=20, s=1.0)
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
kout = krylov_iteration(sn_data,s,1.e-12)
kdiff = norm(ndaout.sol - kout.sol,Inf)
println(kdiff)
semilogy(kout.reshistg)
end
