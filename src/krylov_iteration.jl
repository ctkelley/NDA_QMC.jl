"""
krlov_iteration(sn_data,s,tol=1.e-8; onlygmres=false)
Solves the transport equation with GMRES and BiCGSTAB
"""
function krylov_iteration(sn_data,s,tol=1.e-8; onlygmres=false)
    #
    # Set up the linear system for the Krylov solvers
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
    # kl_gmres and kl_bicgstab solve the problem.
    #
    gout = kl_gmres(f0, rhs, sn_matvec, V, tol; pdata = sn_data)
    if onlygmres
    bout=(solb=[], reshist=[])
    else
    bout = kl_bicgstab(f0, rhs, sn_matvec, V, tol; pdata = sn_data)
    end
    #
    sol = gout.sol
    solb = bout.sol
    reshistg=gout.reshist
    reshistb=bout.reshist
    return(sol=sol, solb=solb, reshistg=reshistg, reshistb=reshistb)
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


