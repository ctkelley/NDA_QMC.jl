function sanity(Nx=2048, na=40, s=1.0)
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
return (snndaok && snok && ndaok)
end
