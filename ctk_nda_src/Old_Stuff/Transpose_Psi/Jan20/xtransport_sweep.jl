"""
transport_sweep!(psi, phi, psi_left, psi_right, angles, source, 
              sweep_data)

Take a single transport sweep.
"""
function xtransport_sweep!(psi, phi, psi_left, psi_right, angles, source, 
              sweep_data)

c=sweep_data.c;
q=source;
dx=sweep_data.dx;
ptmp=sweep_data.ptmp;
na2=length(angles);
na=floor(Int,na2/2);
nx=length(phi);
#
# Initialize the intensities.
#
psi *= 0.0
@views psi[nx,1:na]=psi_right;
@views psi[1,na+1:na2]=psi_left;
#
#
#
source_total = .5*c.*phi + q;
@views source_average=(source_total[2:nx]+source_total[1:nx-1])*.5;
@views forward_angles=angles[na+1:na2];
@views backward_angles=angles[1:na];
#
vfl=(forward_angles/dx) .+ .5;
vfl=1.0 ./ vfl;
vfr=(forward_angles/dx) .- .5;
dfr=Diagonal(vfr);
#
# Forward sweep
#
ptmp.=psi[1,na+1:na2]
@views for ix=2:nx
    psi[ix,na+1:na2] .= ptmp.*vfr
#    psi[ix,na+1:na2] .= psi[ix-1,na+1:na2]
#    psi[ix,na+1:na2] .*= vfr
    psi[ix,na+1:na2] .+= source_average[ix-1]
    psi[ix,na+1:na2] .*= vfl
    ptmp.=psi[ix,na+1:na2]
end
#
#
# Backward sweep
#
ptmp.=psi[nx,1:na]
@views for ix=nx-1:-1:1
    psi[ix,1:na] .= ptmp.*vfr
    psi[ix,1:na] .+= source_average[ix]
    psi[ix,1:na] .*= vfl
    ptmp.=psi[ix,1:na]
end
return psi
end
