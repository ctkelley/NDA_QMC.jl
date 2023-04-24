"""
ctof(uc,uf)

coarse to fine intergrid transfer. This is a cell average deal and
the number of cells must be a power of 2
"""
function ctof(uc,uf)
nf=length(uf)
nc=length(uc)
lnf=Int(log2(nf))
(2^lnf == nf) || error("number of cells must be power of 2")
lnc=Int(log2(nc))
(2^lnc == nc) || error("number of cells must be power of 2")
ncells = Int(nf/nc)
for ic=1:nc
  cbase=(ic - 1) * ncells + 1
  uf[cbase:cbase+ncells-1] .= uc[ic]
end
return uf
end


