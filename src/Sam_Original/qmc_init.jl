using Sobol

function qmc_init(N, Nx, na2, s)

    Lx = 5.0
    dx = Lx/Nx
    #define tally mesh
    low_edges = range(0, stop=Lx-dx, length=Nx)
    high_edges = low_edges.+dx
    dxs = high_edges - low_edges
    midpoints = 0.5*(high_edges + low_edges)
    edges = range(0, stop=Lx, length=Nx+1)

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

    # multi group/zone matrix
    sigt = [1] # pseudo 4-group setup
    G = size(sigt)[1] # number of groups
    temp = zeros(Nx, G)
    for i in 1:G
        temp[:,i] = ones(Nx)*sigt[i]
    end
    sigt = temp
    # source
    source_strength = 0.0
    source = source_strength*ones(Nx,G)
    # scattering
    sigsFunc(x) = exp.(-x/s)
    # s now needs to be input as an array that matches the dimensions of sigt
    sigs = sigsFunc(midpoints)
    #sigs = ones(Nx,G)
    #phi_avg is defaulted to = zeros(Nx)
    phi_edge = zeros(Nx+1,G)
    dphi = zeros(Nx,G)
    phi_s = zeros(Nx,G) .+ 1e-6
    # current
    J_avg = zeros(Nx,G)
    J_edge = zeros(Nx+1,G)

    return qmc_data = (
        N = N,
        Nx = Nx,
        Lx = Lx,
        dx = dx,
        low_edges = low_edges,
        high_edges = high_edges,
        dxs = dxs,
        midpoints = midpoints,
        edges = edges,
        dmu = dmu,
        exit_left_bins = exit_left_bins,
        exit_right_bins = exit_right_bins,
        phi_edge = phi_edge,
        dphi = dphi,
        phi_s = phi_s,
        J_avg = J_avg,
        J_edge = J_edge,
        sigt = sigt,
        source = source,
        sigs = sigs,
        c = s,
        G = G)
end
