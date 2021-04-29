"""
qmc_gmres(N=10^4, Nx=100, na2=11; s=1.0)
"""
function qmc_gmres(N=10^3, Nx=100, na2=11; s=1.0, tol=1.e-10)
mdata=mdata_init(N, Nx, na2, s)
phi0=zeros(Nx,);
b=mdata.frhs;
V=zeros(Nx,20)
eta=tol
gmout=kl_gmres(phi0,b,axb,V,eta; pdata=mdata);
return gmout
end


function tab_test(N=10^3, Nx=100, na2=11; s=1.0, tol=1.e-8)
itout=qmc_gmres(N, Nx, na2; s=s, tol=tol)
qmctab=sn_tabulate(s, Nx, itout.sol; maketab=false, phiedge=false)
return [qmctab.left qmctab.right]
end


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


function axb(phi,mdata)
qmc_data=mdata.qmc_data;
frhs=mdata.frhs
matvec=Qlin(phi,qmc_data,frhs)
end

function Qlin(phi,qmc_data,frhs)
nullout=qmc_sweep(phi,qmc_data)
Mprod=phi-(nullout.phi_avg-frhs)
end

function mdata_init(N, Nx, na2, s)
qmc_data = qmc_init(N, Nx, na2, s);
phizero=zeros(Nx,);
nullout=qmc_sweep(phizero,qmc_data);
frhs=nullout.phi_avg;
mdata=(frhs=frhs, qmc_data =qmc_data)
end
