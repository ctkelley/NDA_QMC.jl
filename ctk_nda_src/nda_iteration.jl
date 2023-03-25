function nda_compare(s=1.0, nx=2001, na=20)
rpic=nda_iteration(nx,na,s);
ndakout=nda_nsoli(nx,na,s);
pich=rpic.reshist./rpic.reshist[1];
lpic=length(pich);
apic=collect(1:lpic);
rnew=ndakout.nout
newh=(rnew.history/rnew.history[1])*pich[2];
newabs= rnew.stats.ifun+ rnew.stats.ijac;
anew=cumsum(newabs);
anew[1]+=1;
bnew=ndakout.noutb
newhb=(bnew.history/bnew.history[1])*pich[2];
newabsb = bnew.stats.ifun+ bnew.stats.ijac;
anewb=cumsum(newabsb);
anewb[1]+=1;
semilogy(apic,pich,"k-",anew,newh,"k--",anewb, newhb,"k-.");
legend(["Picard", "Newton-GMRES", "Newton-Bi-CGSTAB"])
xlabel("Transport Sweeps")
ylabel("Relative Residual")
return (newh=newh, pich=pich, newhb=newhb)
end

function nda_iteration(nx, na, s=1.0)
nda_data=nda_init(nx,na,s)
phi=zeros(nx)
delphi=2.0
itc=1
maxit=200
ithist=zeros(maxit)
while (itc < maxit+1) && (delphi > 1.e-8)
    phiout=nda_fixed(phi,nda_data)
    delphi=norm(phi-phiout,Inf)
    ithist[itc]=delphi
    phi=phiout
    itc += 1
end
#sn_tabulate(s, nx, phi)
reshist=ithist[1:itc-1]
return (reshist=reshist, sol=phi)
end

function nda_nsoli(nx, na, s=1.0)
FS=zeros(nx,)
FPS=zeros(nx,10)
phi0=zeros(nx,)
nda_data=nda_init(nx,na,s)
#
# Fix the initial iterate if you want decent results.
#
phi0=nda_fixed(phi0,nda_data)
nout=nsoli(Fnda!, phi0, FS, FPS; eta=.1, fixedeta=false,
              rtol=1.e-10, pdata=nda_data)
noutb=nsoli(Fnda!, phi0, FS, FPS; eta=.1, fixedeta=false,
              rtol=1.e-10, lsolver="bicgstab", pdata=nda_data)
return (nout=nout, noutb=noutb)
end
    
