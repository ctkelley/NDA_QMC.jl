

function aniso(rn)
    siga = 1
    q0 = 0.75
    q1 = 0.5
    z = 0

    denominator = q0/(2*siga) + z*q1/(2*siga) - q1/(3*siga^2)
    f(mu) = (mu^2*q0/(2*siga) + mu^2*z*q1/(2*siga) - mu^3*q1/(3*siga^2))/denominator - rn

    a = 0
    b = 1
    c = 0
    N=0
    maxiter = 100
    tol = 1e-8
    # bisection method
    while (N<=maxiter) && ((b-a)*0.5 > tol)
        c = (a+b)*0.5
        if ((f(c)<0) && (f(a)<0)) || ((f(c)>0) && (f(a)>0))
            a = c
        else
            b = c
        end
        N += 1
    end
    return c
end

function iso(rn)
    return sqrt(rn)
end

function cdf()
    siga = 1
    q0 = 1.0
    q1 = 0.75
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
