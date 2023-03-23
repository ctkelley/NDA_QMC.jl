"""
sn_angles(na2=40)

Get double Gauss nodes and weights for SN
This function uses FastGaussQuadrature
"""
function sn_angles(na2 = 40)
    na = floor(Int, na2 / 2)
    2 * na == na2 || error("odd number of angles")
    baseangles, baseweights = gausslegendre(na)
    posweights = baseweights * 0.5
    negweights = copy(posweights)
    posangles = (baseangles .+ 1.0) * 0.5
    negangles = -copy(posangles)
    weights = [negweights; posweights]
    angles = [negangles; posangles]
    angles, weights
end
