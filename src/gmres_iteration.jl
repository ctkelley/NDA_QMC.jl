function gmres_iteration(sn_data,s,tol=1.e-8)
    #
    # Set up the linear system for GMRES
    #
    nx=sn_data.nx
    #
    # Separate the right side from the linear operator
    # First step is apply the fixed point map with the real
    # boundary conditions and zero flux.
    #
    rhs = zeros(nx,)
    rhs = flux_map!(rhs, sn_data)
    #
    # Now set the boundary conditions to zero in the
    # precomputed data. That sets up the linear operator.
    #
    sn_data.psi_left .*= 0.0
    sn_data.psi_right .*= 0.0
    #
    # Get organized to call GMRES. Allocate storage for the solution
    # and the Krylov basis
    #
    f0 = zeros(size(rhs))
    V = zeros(nx, 20)
    #
    # kl_gmres solves the problem.
    #
    gout = kl_gmres(f0, rhs, sn_matvec, V, tol; pdata = sn_data)
    #
    sol = gout.sol
    reshist=gout.reshist
    #
    # Tabulate the exit distributions to check results.
    #
    tout = sn_tabulate(s, nx, sol)
    return (left = tout.left, right = tout.right, sol=sol , reshist=reshist)
end

"""
sn_matvec(f, sn_data)

Take the source iteration map and extract the
linear operator part. This is the mat-vec for GMRES.
"""
function sn_matvec(f, sn_data)
    kf = copy(f)
    kf = flux_map!(kf, sn_data)
    kf = f - kf
    return kf
end


