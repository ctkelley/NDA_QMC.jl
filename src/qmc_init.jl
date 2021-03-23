

function qmc_init(Nx, na2, s)

    Lx = 5.0
    dx = Lx/Nx
    #define tally mesh
    low_edges = range(0, stop=Lx-dx, length=Nx)
    high_edges = low_edges.+dx
    dxs = high_edges - low_edges
    midpoints = 0.5*(high_edges + low_edges)
    edges = range(0, stop=Lx, length=Nx+1)

    fl(mu) =  sqrt(mu)
    fr(mu) = -sqrt(mu)

    #initialize sobol sequence
    rng = SobolSeq(2)

    if (N >0)
        skip(rng,N) #skipping the expected number is suggested for Sobol
    end

    rng_bndl = 0
    if (has_left > 0)
        rng_bndl = SobolSeq(2)
        skip(rng_bndl,Nb) #skipping the expected number is suggested for Sobol
    end

    #define angular flux mesh
    #exit_left_bins data structure to hold the exiting angular flux,
    #the first column has the bin centers and the
    #and the second holds the values.
    #exit_right_bins is the same
    dmu = 1/na2
    #right bins
    exit_right_bins = zeros((na2,2))
    exit_right_bins[:,1] = range(dmu/2, stop = 1- dmu/2, step = dmu)
    exit_right_bins[:,2] .= 0
    #left bins
    exit_left_bins = zeros((na2,2))
    exit_left_bins[:,1] = -1*range(dmu/2, stop = 1- dmu/2, step = dmu)
    exit_left_bins[:,2] .= 0

    #phi_avg is defaulted to = zeros(Nx)
    phi_edge = zeros(Nx+1)
    dphi = zeros(Nx)
    phi_s = zeros(Nx) .+ 1e-6

    J_avg = zeros(Nx)
    J_edge = zeros(Nx+1)

    sigt = ones(Nx)*sigt
    source = source_strength*ones(Nx)
    sigsFunc(x) = exp.(-x/c)
    sigs = sigsFunc(midpoints)
    q = phi_avg.*sigs .+ source

    return qmc_data= ()
