"""
sn_init(nx, na2, s; siewert=false, phiedge=true)

I pass a named tuple of precomputed and preallocated data to
all the functions and solvers. 

The input to this is obvious stuff.

nx = number of spatial grid points

na2 = number of angles. The angular mesh is (na2/2) Gaussian quadaratures
      on [-1,0) and (0,1]

s = parameter in Garcia/Siewert

angles, weights = nodes and weights for the angular quadrature
"""
function sn_init(nx, na2, s; siewert=false, phiedge=true)
    if siewert
    angles = [-.05; collect(-.1:-.1:-1.0); 0.05; collect(0.1:0.1:1.0)]
    weights = angles
    else
    (angles, weights) = sn_angles(na2)
    end
    tau = 5.0
    dx = tau / (nx - 1)
    na = floor(Int, na2 / 2)
    x = collect(0:dx:tau)
    c = exp.(-x / s)
    psi_right = zeros(na)
    psi_left = ones(na)
    source = zeros(nx)
    ptmp = zeros(na)
    phiedge ? np=nx : np=nx+1
    psi = zeros(na2, np)
    return sn_data = (
        c = c,
        dx = dx,
        psi = psi,
        angles = angles,
        weights = weights,
        nx = nx,
        source = source,
        psi_right = psi_right,
        psi_left = psi_left,
        ptmp = ptmp,
    )
end
