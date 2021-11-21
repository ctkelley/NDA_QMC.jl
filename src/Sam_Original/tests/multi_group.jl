"""
Sam Pasmann
"""

function qmc_multiGroup()
    ###############################################################################
    #### Parameters
    ###############################################################################

    Nx = 50     # number of tally cells
    na2 = 11    # number of angles for angular mesh
    s = [1]     # parameter in Garcia/Siewert
    N = 2^11    # number of particles per itertion per source
    LB = 0      # left bound
    RB = 500    # right bound
    G = 12      # number of groups
    geometry = "Slab"
    generator = "Sobol"

    qmc_data = multiGroup_init(geometry, generator, N, LB, RB, Nx, G)
    phi_avg, phi_edge, dphi, J_avg, J_edge, history, itt = qmc_source_iteration(s,qmc_data)

    ###############################################################################
    #### Plots
    ###############################################################################
    midpoints = qmc_data.midpoints
    edges = qmc_data.edges
    sol = qmc_data.true_flux

    figure()
    title("Scalar Flux")
    plot(midpoints, sum(phi_avg, dims=2), label=generator)
    plot(edges, sum(phi_edge, dims=2), label="edges")
    ylabel("cell averaged flux")
    xlabel("midpoints")
    legend()


    counter = 1
    figure()
    color_sequence = ["#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
                      "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
                      "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
                      "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"]
    title("Group Scalar Flux")
    for i in 1:1:G
        plot(midpoints, phi_avg[:,i], label=i,color=color_sequence[i] ) #color=color_sequence[i]
        plot(midpoints,sol[i]*ones(Nx),"--",color=color_sequence[i] )
    end
    ylabel("cell averaged flux")
    xlabel("midpoints")
    legend()

end
