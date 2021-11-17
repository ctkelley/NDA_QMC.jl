
using PyPlot
pygui(true)

function bisection(f, a, b, tol=1e-5, maxiter=100)
    fa = f(a)
    fa*f(b) <= 0 || error("No real root in [a,b]")
    i = 0
    local c
    while b-a > tol
        i += 1
        i != maxiter || error("Max iteration exceeded")
        c = (a+b)/2
        fc = f(c)
        if fc == 0
            break
        elseif fa*fc > 0
            a = c  # Root is in the right half of [a,b].
            fa = fc
        else
            b = c  # Root is in the left half of [a,b].
        end
    end
    return c
end


function aniso(rn)
    """
    anisotropic solve.

    given a random number [0,1] solve for the respective mu using
    direct inversion sampling
    """
    siga = 1.0
    q0 = 0.75
    q1 = 0.0
    z = 0

    denominator = q0/(2*siga) + z*q1/(2*siga) - q1/(3*siga^2)
    f(mu) = (mu^2*q0/(2*siga) + mu^2*z*q1/(2*siga) - mu^3*q1/(3*siga^2))/denominator - rn

    a = 0
    b = 1
    c = 0
    N=0
    maxiter = 100
    tol = 1e-6

    mu = bisection(f, a, b, tol, maxiter)

    return mu
end

function iso(rn)
    """
    direct inversion sampling of isotropic source
    """
    return sqrt(rn)
end

function cdf()
    """
    CDF of anisotropic distribution
    """
    siga = 1
    q0 = 1.0
    q1 = 0.0
    z = 0

    denominator = q0/(2*siga) + z*q1/(2*siga) - q1/(3*siga^2)
    f(mu) = (mu.^2*q0./(2*siga) .+ mu.^2*z.*q1./(2*siga) .- mu.^3*q1./(3*siga^2))./denominator
    mus = 0:0.01:1
    figure()
    plot(mus, f(mus))
    title("Anisotropic CDF")

    return
end


cdf()
N = 100000
mus1 = zeros(N)
mus2 = zeros(N)
for i in 1:N
    rn = rand()
    mus1[i] = aniso(rn)
    mus2[i] = iso(rn)
end

bins = 20
figure()
subplot(122)
hist(mus1,bins=bins)
title("anisotropic")
subplot(121)
hist(mus2, bins=bins)
title("isotropic")
