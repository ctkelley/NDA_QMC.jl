"""
sn_angles(na=20)

Get double Gauss nodes and weights for SN
"""
function sn_angles(na=20)
na2=floor(Int,na/2)
2*na2 == na || error("odd number of angles")
baseangles, baseweights = gausslegendre(na2)
posweights=baseweights * .5
negweights=copy(posweights)
posangles=(baseangles .+ 1.0) *.5
negangles=-copy(posangles)
weights=[negweights; posweights]
angles=[negangles; posangles]
angles, weights
end

