"""
flux_map
"""
function flux_map!(flux,psi_left,psi_right,sn_data)
psi=sn_data.psi
weights=sn_data.weights
psi = transport_sweep!(psi, flux, psi_left, psi_right, sn_data);
flux.= (weights'*psi)'
end

