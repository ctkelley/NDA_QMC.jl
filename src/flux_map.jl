"""
flux_map
"""
function flux_map!(flux,sn_data)
psi=sn_data.psi
weights=sn_data.weights
psi_left=sn_data.psi_left
psi_right=sn_data.psi_right
psi = transport_sweep!(psi, flux, sn_data);
flux.= (weights'*psi)'
end

