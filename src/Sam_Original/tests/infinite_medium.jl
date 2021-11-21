
function qmc_infMed()
    ###############################################################################
    #### Parameters
    ###############################################################################

    Nx = 40     # number of tally cells
    N = 2^11    # number of particles per source itertion
    LB = 0      # left bound
    RB = 10     # right bound
    sigt = [1.0]
    sigs = 0.5*sigt
    siga = sigt - sigs
    geometry = "Slab"   # infinite medium problems are only configured for slab geometries right now
    generator = "Sobol"

    ###############################################################################
    #### Function Call
    ###############################################################################

    qmc_data = const_infMed_init(geometry, generator, N, LB, RB, Nx, sigs, sigt)
    @time begin
    phi_avg, phi_edge, dphi, J_avg, J_edge, psi_right, psi_left, history, itt = qmc_source_iteration(qmc_data,1.e-3)
    end

    ###############################################################################
    #### Plotting
    ###############################################################################
    midpoints = qmc_data.midpoints
    flux = qmc_data.phi_avg
    sol = qmc_data.phi_true

    edges = qmc_data.edges

    figure()
    plot(midpoints, phi_avg, label = generator)
    plot(edges, phi_edge, label="edges")
    plot(midpoints, sol, label= "solution")
    ylabel("Cell Average Flux")
    xlabel("Cell Midpoints")
    title("Scalar Flux")
    legend()

end
