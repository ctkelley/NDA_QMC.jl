"""
flux_map

This is the source iteration fixed point map. The input is

flux = right side of the fixed source problem
sn_data = the precomputed data for this

Notice that I've preallocated the storage for the angular flux psi.
I've stored it as an na2 by nx array so the flux computation has
a few transposes in it.

The output is the scalar flux (and current when the times comes)
you get from the fixed source solve.

"""
function flux_map!(flux, sn_data)
    psi = sn_data.psi
    psi = transport_sweep!(psi, flux, sn_data)
    #
    # Take the 0th moment to get the flux.
    #
    weights = sn_data.weights
    flux .= (weights' * psi)'
end
