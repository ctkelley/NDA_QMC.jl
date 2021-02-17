"""
sn_tabulate(s,nx,flux)

Make the tables to compare with Garcia/Siewert

Uses the converged flux from the solve.
"""
function sn_tabulate(s, nx, flux)
    angleout = [-.05; collect(-.1:-.1:-1.0); 0.05; collect(0.1:0.1:1.0)]
    #
    # I don't really need the weights, but sn_init expects some
    weights = angleout
    #
    na2 = length(angleout)
    na = floor(Int, na2 / 2)
#    sn_data = sn_init(nx, na2, s, angleout, weights)
    sn_data = sn_init(nx, na2, s; siewert=true)
    psi = sn_data.psi
    psi = transport_sweep!(psi, flux, sn_data)
    header = " mu         I(0,-mu)        I(tau,mu)"
    @printf("%s \n", header)
    for it=1:na
    @printf("%5.2f %15.5e %15.5e \n", 
           angleout[it+na], psi[it,1], psi[na+it,nx])
    end
    return (left = psi[1:na, 1], right = psi[na+1:na2, nx])
end
