"""
qmc_krylov(N=10^4, Nx=100, na2=11; s=1.0)
Solves the linear system you get from QMC with GMRES.
"""
function qmc_krylov(N=10^3, Nx=100, na2=11; s=1.0, tol=1.e-10, onlygmres=true)
mdata=mdata_init(N, Nx, na2, s)
phi0=zeros(Nx,);
b=mdata.frhs;
V=zeros(Nx,20)
eta=tol
#
# kl_gmres and kl_bicgstab solve the problem.
#
gout = kl_gmres(phi0, b, axb, V, tol; pdata = mdata)
if onlygmres
    bout=(sol=[], reshist=[])
    else
    bout = kl_bicgstab(phi0, b, axb, V, tol; pdata = mdata)
end
#gmout=kl_gmres(phi0,b,axb,V,eta; pdata=mdata);
sol=gout.sol
solb=bout.sol
reshistg=gout.reshist
reshistb=bout.reshist
return(sol=sol, solb=solb, reshistg=reshistg, reshistb=reshistb)
end


"""
tab_test(N=10^3, Nx=100, na2=11; s=1.0, tol=1.e-8)
Makes a table to compare to Garcia-Siewert
"""
function tab_test(N=10^3, Nx=100, na2=11; s=1.0, tol=1.e-8)
itout=qmc_krylov(N, Nx, na2; s=s, tol=tol, onlygmres=true)
qmctab=sn_tabulate(s, Nx, itout.sol; maketab=false, phiedge=false)
return [qmctab.left qmctab.right]
end


"""
qmc_si(N=10^4, Nx=100, na2=11; s=1.0, tol=1.e-8, maxit=100)
Source iteration using QMC. Nothing magic here.
"""
function qmc_si(N=10^4, Nx=100, na2=11; s=1.0, tol=1.e-8, maxit=100)
qmc_data=qmc_init(N, Nx, na2, s)
phic=zeros(Nx,);
delflux=1.0
reshist=[]
itt=1
while itt < maxit && delflux > tol
nullout=qmc_sweep(phic,qmc_data)
phip = nullout.phi_avg
delflux=norm(phip-phic,2)
itt += 1
push!(reshist, delflux)
phic .= phip
end
return ( sol=phic, history=reshist)
end


"""
axb(phi,mdata)
Matrix-vector product to feed to linear solvers.
"""
function axb(phi,mdata)
qmc_data=mdata.qmc_data;
frhs=mdata.frhs
matvec=Qlin(phi,qmc_data,frhs)
end

"""
Qlin(phi,qmc_data,frhs)
Does a QMC transport sweep with flux phi, subtracts off the zero bc solution
to obtain the linear integral operator. The subtract that from the identity. 
"""
function Qlin(phi,qmc_data,frhs)
nullout=qmc_sweep(phi,qmc_data)
Mprod=phi-(nullout.phi_avg-frhs)
end


"""
mdata_init(N, Nx, na2, s)
Collects the precomputed data for QMC and does a sweep with 
zero boundary data to get the right hand side for the linear
equation forumlation.
"""
function mdata_init(N, Nx, na2, s)
# Precomputed data
qmc_data = qmc_init(N, Nx, na2, s);
# Sweep with zero RHS
phizero=zeros(Nx,);
nullout=qmc_sweep(phizero,qmc_data);
frhs=nullout.phi_avg;
#
mdata=(frhs=frhs, qmc_data =qmc_data)
end
