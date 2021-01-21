"""
sn_tabulate(s,nx,flux)

Make the tables to compare with Garcia/Siewert

Uses the converged flux from the solve.
"""
function sn_tabulate(s,nx,flux)
angleout=[-.05; collect(-.1:-.1:-1.0); .05; collect(.1:.1:1.0)]
weights=angleout
na2=length(angleout)
na=floor(Int,na2/2)
sn_data=sn_init(nx, na2, s, angleout, weights)
psi=sn_data.psi
psi = transport_sweep!(psi, flux, sn_data);
return (left=psi[1:na,1], right=psi[na+1:na2,nx])
end



