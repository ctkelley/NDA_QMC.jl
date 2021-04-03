"""
siewert(s=1.0)

Gets the Garcia-Siewert tabular data and parks it where it's easy to find.
"""
function siewert(s=1.0; filedump=false)
    na2 = 20
    nx = 2001
    #
    # precomputed data
    #
    sn_data = sn_init(nx, na2, s)
    #
    # Solve the problem
    #
    gout=gmres_iteration(sn_data,s)
    fgm = gout.sol
    tabout=sn_tabulate(s, nx, fgm; maketab=false)
    Dout=[tabout.left tabout.right]
    if filedump
    s == 1.0 ? sstr=string(1) : sstr="Inf"
    fname=string("Siewert:s=",sstr,".dat")
    writetab(fname,Dout)
    end
    return Dout
end    

function readtab(fname,data)
open(fname,"r") do io
read!(io,data)
end
end

function writetab(fname,data)
open(fname,"w") do io
write(io, data)
end
end
