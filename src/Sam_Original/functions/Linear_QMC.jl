"""
Linear_QMC()

Demonstrate a lack of linearity. I start with the issue in QMC_Sweep
and show that the obvious way to recover linear maps does not work.

See the pdf for the details.

"""

function Linear_QMC()
    qmc_data = EasyInit()
    b = getb(qmc_data)
    Nx = qmc_data.Nx
    phi0 = ones(Nx)
    phi1 = TSweep(phi0, qmc_data)
    r0 = phi1 - phi0
    #
    # Linearity?
    #
    v1 = phi0 + r0
    mv1 = MMul(v1, qmc_data, b)
    mphi0 = MMul(phi0, qmc_data, b)
    mr0 = MMul(r0, qmc_data, b) # this is zero!
    linerror = mv1 - (mphi0 + mr0) # supposed to be zero
    println("This ", norm(linerror), " is supposed to be zero.")
end

function test2()
    qmc_data = EasyInit()
    b = getb(qmc_data)
    Nx = qmc_data.Nx
    phi0 = ones(Nx)
    phi1 = TSweep(phi0, qmc_data)
    r0 = phi1 - phi0
    mr0 = MMul(r0, qmc_data, b)
    # mr0=0. That's the problem I sent yesterday.
    println("Mr0 =", mr0, " Here's the norm ", norm(mr0))
end


function MMul(phi, qmc_data, b)
    mmul = TSweep(phi, qmc_data) - b
    return mmul
end

function test1()
    qmc_data = EasyInit()
    b = getb(qmc_data)
    Nx = qmc_data.Nx
    linerror = 0.0
    for itest = 1:100
        x = rand(Nx)
        y = rand(Nx)
        alpha = rand()
        beta = rand()
        z = alpha * x + beta * y
        mz = MMul(z, qmc_data, b)
        mx = MMul(x, qmc_data, b)
        my = MMul(y, qmc_data, b)
        mxy = alpha * mx + beta * my
        linerror += norm(mz - mxy)
    end
    return linerror
end


function TSweep(phi, qmc_data)
    phiout = qmc_sweep(phi, qmc_data)
    phinew = phiout.phi_avg
    return phinew
end

function getb(qmc_data)
    Nx = qmc_data.Nx
    zed = zeros(Nx)
    b = TSweep(zed, qmc_data)
    return b
end


function EasyInit()
    #
    # Init with usual stuff
    #
    Nx = 10     # number of tally cells
    na2 = 11    # number of angles for angular mesh
    s = 1.0     # parameter in Garcia/Siewert
    N = 100   # number of particles per source itertion
    LB = 0.0   # left bound
    RB = 5.0   # right bound
    geometry = "Slab"
    generator = "Sobol"
    qmc_data = garcia_init(geometry, generator, N, LB, RB, Nx, na2, s)
end
