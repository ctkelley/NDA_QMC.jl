"""
transport_sweep!(psi, phi, psi_left, psi_right, sn_data)

Take a single transport sweep.
"""
function transport_sweep!(psi, phi, psi_left, psi_right, sn_data)

angles=sn_data.angles
method = sn_data.method
source = sn_data.source
#
ptmp=sn_data.ptmp
c=sn_data.c;
dx=sn_data.dx;
#
na2=length(angles);
na=floor(Int,na2/2);
nx=length(phi);
#
# Initialize the intensities.
#
psi *= 0.0
if method == :old
@views psi[nx,1:na]=psi_right;
@views psi[1,na+1:na2]=psi_left;
else
@views psi[1:na,nx]=psi_right;
@views psi[na+1:na2,1]=psi_left;
end
#
#
#
source_total = .5*c.*phi + source;
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
if method == :old
ptmp.=psi[1,na+1:na2]
@views for ix=2:nx
    psi[ix,na+1:na2] .= ptmp.*vfr
    psi[ix,na+1:na2] .+= source_average[ix-1]
    psi[ix,na+1:na2] .*= vfl
    ptmp.=psi[ix,na+1:na2]
end
else
#ptmp.=psi[na+1:na2,1]
@views for ix=2:nx
    copy!(psi[na+1:na2,ix],psi[na+1:na2,ix-1])
    psi[na+1:na2,ix] .*= vfr
    psi[na+1:na2,ix] .+= source_average[ix-1]
    psi[na+1:na2,ix] .*= vfl
end

end
#
#
# Backward sweep
#
if method == :old
ptmp.=psi[nx,1:na]
@views for ix=nx-1:-1:1
    psi[ix,1:na] .= ptmp.*vfr
    psi[ix,1:na] .+= source_average[ix]
    psi[ix,1:na] .*= vfl
    ptmp.=psi[ix,1:na]
end
else
#ptmp.=psi[1:na,nx]
@views for ix=nx-1:-1:1
    copy!(psi[1:na,ix],psi[1:na,ix+1])
    psi[1:na,ix] .*= vfr
    psi[1:na,ix] .+= source_average[ix]
    psi[1:na,ix] .*= vfl
end
end
return psi
end
