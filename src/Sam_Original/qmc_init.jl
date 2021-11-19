using DelimitedFiles
using LinearAlgebra
"""
The functions below used to precompute data for various problem types.
"""

"""
    qmc_init(Geometry, generator, N, LB, RB, Nx, na2, s, sigs, sigt)
General init function not specific to any problem type. Originally used for Garcia/
Siewart problems.
"""
function qmc_init(Geometry, generator, N, LB, RB, Nx, na2, s, sigs, sigt)

        hasLeft = false
        hasRight = false

        dx = (RB-LB)/Nx
        #define tally mesh
        low_edges = range(LB, stop=RB-dx, length=Nx)
        high_edges = low_edges.+dx
        midpoints = 0.5*(high_edges + low_edges)
        edges = range(LB, stop=RB, length=Nx+1)

        #define angular flux mesh
        #exit_left_bins data structure to hold the exiting angular flux,
        #the first column has the bin centers and
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

        G = size(sigt)[1] # number of groups
        temp1 = zeros(Nx, G)
        for i in 1:G
            temp1[:,i] = ones(Nx)*sigt[i]
        end
        sigt = temp1
        source_strength = 1.0
        source = source_strength*ones(Nx,G)

        # phi_avg is defaulted to = zeros(Nx)
        phi_edge = zeros(Nx+1,G)
        phi_avg = source_strength*zeros(Nx,G)
        # used for boundary sources
        phi_left = 1.0  # source_strength
        phi_right = 1.0 # source_strength

        dphi = zeros(Nx,G)
        phi_s = zeros(Nx,G) .+ 1e-6
        # current
        J_avg = zeros(Nx,G)
        J_edge = zeros(Nx+1,G)
        # Garcia parameter
        c = s
        # Geometry
        if (Geometry == "Slab")
            Geo = 1
        elseif (Geometry == "Cylinder")
            Geo = 2
        elseif (Geometry == "Sphere")
            Geo = 3
        else
            error("Geometry must be: 'Slab', 'Cylinder', or 'Sphere'")
        end

        qmc_data = (        Geo = Geo,
                            N = N,
                            Nx = Nx,
                            G = G,
                            RB = RB,
                            LB = LB,
                            dmu = dmu,
                            exit_left_bins = exit_left_bins,
                            exit_right_bins = exit_right_bins,
                            phi_edge = phi_edge,
                            phi_avg = phi_avg,
                            phi_left = phi_left,
                            phi_right = phi_right,
                            dphi = dphi,
                            phi_s = phi_s,
                            J_avg = J_avg,
                            J_edge = J_edge,
                            sigt = sigt,
                            sigs = sigs,
                            source = source,
                            c = c,
                            generator = generator,
                            hasRight = hasRight,
                            hasLeft = hasLeft,
                            low_edges = low_edges,
                            high_edges = high_edges,
                            midpoints = midpoints,
                            edges = edges)
        return qmc_data
end

"""
    multiGroup_init(Geometry, generator, N, LB, RB, Nx, numGroups)
Multi-Group initialization function. Solution data stored in ./HDPE_data for 12,
70, and 618 groups.
"""
function multiGroup_init(Geometry, generator, N, LB, RB, Nx, numGroups)

        if (numGroups != 12) && (numGroups != 70) && (numGroups != 618)
            error("Number of groups must be 12, 70, or 618.")
        end
        D = readdlm(joinpath(@__DIR__, "./HDPE_data/D_$(numGroups)G_HDPE.csv"), ',', Float64)
        siga = readdlm(joinpath(@__DIR__, "./HDPE_data/Siga_$(numGroups)G_HDPE.csv"), ',', Float64)
        sigs = readdlm(joinpath(@__DIR__, "./HDPE_data/Scat_$(numGroups)G_HDPE.csv"), ',', Float64)
        sigs = reverse(sigs, dims=2)
        sigt = 1 ./ (3*D)

        hasLeft = false
        hasRight = false

        dx = (RB-LB)/Nx
        #define tally mesh
        low_edges = range(LB, stop=RB-dx, length=Nx)
        high_edges = low_edges.+dx
        midpoints = 0.5*(high_edges + low_edges)
        edges = range(LB, stop=RB, length=Nx+1)

        #define angular flux mesh
        dmu = 1/na2
        #right bins
        exit_right_bins = zeros((na2,2))
        exit_right_bins[:,1] = range(dmu/2, stop = 1- dmu/2, step = dmu)
        exit_right_bins[:,2] .= 0
        #left bins
        exit_left_bins = zeros((na2,2))
        exit_left_bins[:,1] = -1*range(dmu/2, stop = 1- dmu/2, step = dmu)
        exit_left_bins[:,2] .= 0

        G = numGroups # number of groups
        source_strength = 1.0
        # analytic solution
        Q = source_strength*ones(numGroups) # source
        true_flux = inv(Diagonal(sigt[:,1]) .- sigs)*Q # returns diagonal matrix

        temp1 = zeros(Nx, G)
        for i in 1:G
            temp1[:,i] = ones(Nx)*sigt[i]
        end
        sigt = temp1
        source = source_strength*ones(Nx,G)

        # phi_avg is defaulted to = zeros(Nx)
        phi_edge = zeros(Nx+1,G)
        phi_avg = source_strength*zeros(Nx,G)
        # used for boundary sources
        phi_left = 1.0  # source_strength
        phi_right = 1.0 # source_strength

        dphi = zeros(Nx,G)
        phi_s = zeros(Nx,G) .+ 1e-6
        # current
        J_avg = zeros(Nx,G)
        J_edge = zeros(Nx+1,G)
        # Garcia parameter
        c = s
        # Geometry
        if (Geometry == "Slab")
            Geo = 1
        elseif (Geometry == "Cylinder")
            Geo = 2
        elseif (Geometry == "Sphere")
            Geo = 3
        else
            error("Geometry must be: 'Slab', 'Cylinder', or 'Sphere'")
        end

        qmc_data = (        Geo = Geo,
                            N = N,
                            Nx = Nx,
                            G = G,
                            RB = RB,
                            LB = LB,
                            dmu = dmu,
                            exit_left_bins = exit_left_bins,
                            exit_right_bins = exit_right_bins,
                            phi_edge = phi_edge,
                            phi_avg = phi_avg,
                            phi_left = phi_left,
                            phi_right = phi_right,
                            dphi = dphi,
                            phi_s = phi_s,
                            J_avg = J_avg,
                            J_edge = J_edge,
                            sigt = sigt,
                            sigs = sigs,
                            source = source,
                            c = c,
                            generator = generator,
                            hasRight = hasRight,
                            hasLeft = hasLeft,
                            low_edges = low_edges,
                            high_edges = high_edges,
                            midpoints = midpoints,
                            edges = edges,
                            true_flux = true_flux)
        return qmc_data
end

"""
    garcia_init(Geometry, generator, N, LB, RB, Nx, na2, s)
Problem set up for Garcia/Siewart paper. Characterized by left boundary source
and spatially dependent scattering cross section.
JQSRT (27), 1982 pp 141-148.
"""
function garcia_init(Geometry, generator, N, LB, RB, Nx, na2, s)

        hasLeft = true
        hasRight = false

        dx = (RB-LB)/Nx
        #define tally mesh
        low_edges = range(LB, stop=RB-dx, length=Nx)
        high_edges = low_edges.+dx
        midpoints = 0.5*(high_edges + low_edges)
        edges = range(LB, stop=RB, length=Nx+1)

        #define angular flux mesh
        #exit_left_bins data structure to hold the exiting angular flux,
        #the first column has the bin centers and
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

        sigt = [1]
        G = size(sigt)[1] # number of groups
        sigt = ones(Nx,G).*sigt
        sigt = hcat(sigt)
        source_strength = 0.0
        source = source_strength*ones(Nx,G)
        sigsFunc(x) = exp.(-x/s)
        sigs = sigsFunc(midpoints)

        #phi_avg is defaulted to = zeros(Nx)
        phi_edge = zeros(Nx+1,G)
        phi_avg = source_strength*zeros(Nx,G)
        dphi = zeros(Nx,G)
        phi_s = zeros(Nx,G) .+ 1e-6
        # current
        J_avg = zeros(Nx,G)
        J_edge = zeros(Nx+1,G)
        # Garcia parameter
        c = s
        # Geometry
        if (Geometry == "Slab")
            Geo = 1
        elseif (Geometry == "Cylinder")
            Geo = 2
        elseif (Geometry == "Sphere")
            Geo = 3
        else
            error("Geometry must be: 'Slab', 'Cylinder', or 'Sphere'")
        end

        phi_left = 1
        phi_right = 0

        qmc_data = (        Geo = Geo,
                            N = N,
                            Nx = Nx,
                            G = G,
                            RB = RB,
                            LB = LB,
                            dmu = dmu,
                            exit_left_bins = exit_left_bins,
                            exit_right_bins = exit_right_bins,
                            phi_edge = phi_edge,
                            phi_avg = phi_avg,
                            phi_left = phi_left,
                            phi_right = phi_right,
                            dphi = dphi,
                            phi_s = phi_s,
                            J_avg = J_avg,
                            J_edge = J_edge,
                            sigt = sigt,
                            sigs = sigs,
                            source = source,
                            c = c,
                            generator = generator,
                            hasRight = hasRight,
                            hasLeft = hasLeft,
                            low_edges = low_edges,
                            high_edges = high_edges,
                            midpoints = midpoints,
                            edges = edges)
        return qmc_data
end

"""
    const_infMed_init(Geometry, generator, N, LB, RB, Nx, na2, s, sigs, sigt)
Infinite medium problem with a constant volumetric source.
"""
function const_infMed_init(Geometry, generator, N, LB, RB, Nx, na2, s, sigs, sigt)

        hasLeft = true
        hasRight = true

        siga = sigt - sigs

        dx = (RB-LB)/Nx
        #define tally mesh
        low_edges = range(LB, stop=RB-dx, length=Nx)
        high_edges = low_edges.+dx
        midpoints = 0.5*(high_edges + low_edges)
        edges = range(LB, stop=RB, length=Nx+1)

        #define angular flux mesh
        #exit_left_bins data structure to hold the exiting angular flux,
        #the first column has the bin centers and
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

        G = size(sigt)[1] # number of groups
        temp1 = zeros(Nx, G)
        for i in 1:G
            temp1[:,i] = ones(Nx)*sigt[i]
        end
        sigt = temp1
        source_strength = 1.0
        source = source_strength*ones(Nx,G)
        phi_true = source./siga
        # phi_avg is defaulted to = zeros(Nx)
        phi_edge = zeros(Nx+1,G)
        phi_avg = source_strength*zeros(Nx,G)

        dphi = zeros(Nx,G)
        phi_s = zeros(Nx,G) .+ 1e-6
        # current
        J_avg = zeros(Nx,G)
        J_edge = zeros(Nx+1,G)
        # Garcia parameter
        c = s
        # Geometry
        if (Geometry == "Slab")
            Geo = 1
        elseif (Geometry == "Cylinder")
            Geo = 2
        elseif (Geometry == "Sphere")
            Geo = 3
        else
            error("Geometry must be: 'Slab', 'Cylinder', or 'Sphere'")
        end

        # used for boundary sources
        phi_left = 0.5*source_strength/siga#*N/surfaceArea(Geo,LB)
        phi_right = 0.5*source_strength/siga#*N/surfaceArea(Geo,RB)

        qmc_data = (        Geo = Geo,
                            N = N,
                            Nx = Nx,
                            G = G,
                            RB = RB,
                            LB = LB,
                            dmu = dmu,
                            exit_left_bins = exit_left_bins,
                            exit_right_bins = exit_right_bins,
                            phi_true = phi_true,
                            phi_edge = phi_edge,
                            phi_avg = phi_avg,
                            phi_left = phi_left,
                            phi_right = phi_right,
                            dphi = dphi,
                            phi_s = phi_s,
                            J_avg = J_avg,
                            J_edge = J_edge,
                            sigt = sigt,
                            sigs = sigs,
                            source = source,
                            c = c,
                            generator = generator,
                            hasRight = hasRight,
                            hasLeft = hasLeft,
                            low_edges = low_edges,
                            high_edges = high_edges,
                            midpoints = midpoints,
                            edges = edges)
        return qmc_data
end

"""
    linear_infMed_init(Geometry, generator, N, LB, RB, Nx, na2, s, sigs, sigt)
Infinite medium problem with a spatially linear dependent source.
"""
function linear_infMed_init(Geometry, generator, N, LB, RB, Nx, na2, s, sigs, sigt)

        hasLeft = true
        hasRight = true

        siga = sigt - sigs

        dx = (RB-LB)/Nx
        #define tally mesh
        low_edges = range(LB, stop=RB-dx, length=Nx)
        high_edges = low_edges.+dx
        midpoints = 0.5*(high_edges + low_edges)
        edges = range(LB, stop=RB, length=Nx+1)

        #define angular flux mesh
        #exit_left_bins data structure to hold the exiting angular flux,
        #the first column has the bin centers and
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

        G = size(sigt)[1] # number of groups
        temp1 = zeros(Nx, G)
        for i in 1:G
            temp1[:,i] = ones(Nx)*sigt[i]
        end
        sigt = temp1

        q0 = 0.75
        q1 = 0.5
        source_strength(x) = q0 .+ q1.*x
        source = source_strength(midpoints).*ones(Nx,G)
        solution(x) = source_strength(x)./siga
        phi_true = solution(midpoints)
        # phi_avg is defaulted to = zeros(Nx)
        phi_edge = zeros(Nx+1,G)
        phi_avg = source_strength(midpoints).*zeros(Nx,G)

        dphi = zeros(Nx,G)
        phi_s = zeros(Nx,G) .+ 1e-6
        # current
        J_avg = zeros(Nx,G)
        J_edge = zeros(Nx+1,G)
        # Garcia parameter
        c = s
        # Geometry
        if (Geometry == "Slab")
            Geo = 1
        elseif (Geometry == "Cylinder")
            Geo = 2
        elseif (Geometry == "Sphere")
            Geo = 3
        else
            error("Geometry must be: 'Slab', 'Cylinder', or 'Sphere'")
        end

        # used for boundary sources
        phi_left = solution(low_edges[1])*0.5
        phi_right = solution(high_edges[Nx])*0.5

        qmc_data = (        Geo = Geo,
                            N = N,
                            Nx = Nx,
                            G = G,
                            RB = RB,
                            LB = LB,
                            dmu = dmu,
                            exit_left_bins = exit_left_bins,
                            exit_right_bins = exit_right_bins,
                            phi_true = phi_true,
                            phi_edge = phi_edge,
                            phi_avg = phi_avg,
                            phi_left = phi_left,
                            phi_right = phi_right,
                            dphi = dphi,
                            phi_s = phi_s,
                            J_avg = J_avg,
                            J_edge = J_edge,
                            sigt = sigt,
                            sigs = sigs,
                            siga = siga,
                            source = source,
                            c = c,
                            generator = generator,
                            hasRight = hasRight,
                            hasLeft = hasLeft,
                            low_edges = low_edges,
                            high_edges = high_edges,
                            midpoints = midpoints,
                            edges = edges)
        return qmc_data
end

"""
    quadratic_infMed_init(Geometry, generator, N, LB, RB, Nx, na2, s, sigs, sigt)
Infinite medium problem with a spatially quadratic source.
"""
function quadratic_infMed_init(Geometry, generator, N, LB, RB, Nx, na2, s, sigs, sigt)

        hasLeft = true
        hasRight = true

        siga = sigt - sigs

        dx = (RB-LB)/Nx
        #define tally mesh
        low_edges = range(LB, stop=RB-dx, length=Nx)
        high_edges = low_edges.+dx
        midpoints = 0.5*(high_edges + low_edges)
        edges = range(LB, stop=RB, length=Nx+1)

        #define angular flux mesh
        #exit_left_bins data structure to hold the exiting angular flux,
        #the first column has the bin centers and
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

        G = size(sigt)[1] # number of groups
        temp1 = zeros(Nx, G)
        for i in 1:G
            temp1[:,i] = ones(Nx)*sigt[i]
        end
        sigt = temp1

        q0 = 0.5
        q1 = 0.75
        q2 = 1.25
        source_strength(x) = q0 .+ q1.*x .+ q2.*x.^2
        source = source_strength(midpoints).*ones(Nx,G)
        solution(x) = (source_strength(x) .+ 2*q2./(3*siga.^2))./siga
        phi_true = solution(midpoints)

        # phi_avg is defaulted to = zeros(Nx)
        phi_edge = zeros(Nx+1,G)
        phi_avg = source_strength(midpoints).*zeros(Nx,G)

        dphi = zeros(Nx,G)
        phi_s = zeros(Nx,G) .+ 1e-6
        # current
        J_avg = zeros(Nx,G)
        J_edge = zeros(Nx+1,G)
        # Garcia parameter
        c = s
        # Geometry
        if (Geometry == "Slab")
            Geo = 1
        elseif (Geometry == "Cylinder")
            Geo = 2
        elseif (Geometry == "Sphere")
            Geo = 3
        else
            error("Geometry must be: 'Slab', 'Cylinder', or 'Sphere'")
        end

        # used for boundary sources
        phi_left = 0.5*solution(low_edges[1])
        phi_right = 0.5*solution(high_edges[Nx])

        qmc_data = (        Geo = Geo,
                            N = N,
                            Nx = Nx,
                            G = G,
                            RB = RB,
                            LB = LB,
                            dmu = dmu,
                            exit_left_bins = exit_left_bins,
                            exit_right_bins = exit_right_bins,
                            phi_true = phi_true,
                            phi_edge = phi_edge,
                            phi_avg = phi_avg,
                            phi_left = phi_left,
                            phi_right = phi_right,
                            dphi = dphi,
                            phi_s = phi_s,
                            J_avg = J_avg,
                            J_edge = J_edge,
                            sigt = sigt,
                            sigs = sigs,
                            source = source,
                            c = c,
                            generator = generator,
                            hasRight = hasRight,
                            hasLeft = hasLeft,
                            low_edges = low_edges,
                            high_edges = high_edges,
                            midpoints = midpoints,
                            edges = edges)
        return qmc_data
end
