"""
transport_sweep!(psi, phi, sn_data; phiedge=true)

Take a single transport sweep.
"""
function transport_sweep!(psi, phi, sn_data; phiedge=true)

    angles = sn_data.angles
    source = sn_data.source
    psi_right = sn_data.psi_right
    psi_left = sn_data.psi_left
    #
    ptmp = sn_data.ptmp
    c = sn_data.c
    dx = sn_data.dx
    #
    na2 = length(angles)
    na = floor(Int, na2 / 2)
    nx = length(phi)
    phiedge ? np=nx : np=nx+1
    #
    # Initialize the angular flux and set the boundary conditions.
    #
    psi *= 0.0
    @views psi[1:na, np] = psi_right
    @views psi[na+1:na2, 1] = psi_left
    #
    #
    #
    source_total = 0.5 * c .* phi + source
    if phiedge
    @views source_average = (source_total[2:nx] + source_total[1:nx-1]) * 0.5
    else
    source_average=source_total
    end
    @views forward_angles = angles[na+1:na2]
    @views backward_angles = angles[1:na]
    #
    vfl = (forward_angles / dx) .+ 0.5
    vfl = 1.0 ./ vfl
    vfr = (forward_angles / dx) .- 0.5
    #
    # Forward sweep
    #
    @views for ix = 2:np
        copy!(psi[na+1:na2, ix], psi[na+1:na2, ix-1])
        psi[na+1:na2, ix] .*= vfr
        psi[na+1:na2, ix] .+= source_average[ix-1]
        psi[na+1:na2, ix] .*= vfl
    end
    #
    # Backward sweep
    #
    @views for ix = np-1:-1:1
        copy!(psi[1:na, ix], psi[1:na, ix+1])
        psi[1:na, ix] .*= vfr
        psi[1:na, ix] .+= source_average[ix]
        psi[1:na, ix] .*= vfl
    end
    return psi
end
