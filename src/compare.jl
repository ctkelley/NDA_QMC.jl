function compare(s=1.0)
na2=20
nx=2001
method=:new
#
# precomputed data
#
(angles, weights) = sn_angles(na2)
sn_data=sn_init(nx,na2,s,angles,weights,method)
sout=source_iteration(na2,s)
fout=sout.flux
rhs=zeros(size(fout))
rhs=flux_map!(rhs,sn_data)
sn_data.psi_left .*=0.0
sn_data.psi_right .*=0.0
kf=copy(fout)
kf=flux_map!(kf,sn_data)
plot(fout-kf-rhs)
f0=zeros(size(fout))
V=zeros(nx,20)
gout=kl_gmres(f0, rhs, sn_matvec, V, 1.e-10; pdata=sn_data);
fgm=gout.sol
figure(1)
plot(fout-fgm)
figure(2)
snhist=sout.history./sout.history[1]
ghist=gout.reshist./gout.reshist[1]
semilogy(ghist,"k-",snhist,"k--")
legend(["gmres", "source iteration"])
xlabel("Transport Sweeps")
ylabel("Relative Residual")
s == Inf ? strs=L"\infty" : strs=string(s)
tstring=string("Residual Histories: s= ",strs)
title(tstring)
println(norm(fout  - fgm,Inf))
end

function sn_matvec(f,sn_data)
kf=copy(f)
kf=flux_map!(kf,sn_data)
kf=f-kf
return kf
end
