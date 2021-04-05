"""
This makes the table in the notebook showing the errors in QMC as a 
function of Nx and N

qmc_vs_sn(levels=6, Nrange=4; s=1.0)
"""
function qmc_vs_sn(levels=6, Nrange=4; s=1.0)
NVals=zeros(Int64,Nrange)
NVals[1]=1000
for i=2:Nrange
NVals[i]=2*NVals[i-1]
end
NxVals=zeros(Int64,levels);
Dcomp=zeros(levels, Nrange);
#
NxVals[1]=50;
for il=2:levels
NxVals[il]=2*NxVals[il-1];
end
DataGS=readdata(s)
for idv=1:Nrange
Dcomp[:,idv]=validate(NVals[idv], levels, DataGS, s)
end
makeqmctab(Dcomp,NVals,NxVals)
return Dcomp
end

function sn_validate(Nrange=3; s=1.0, phiedge=true)
NVals=zeros(Int64,Nrange)
NVals[1]=100
for i=2:Nrange
NVals[i]=2*NVals[i-1]
end
Dcomp=zeros(Nrange);
#
DataGS=readdata(s)
for idv=1:Nrange
Dval=siewert(s; nx=NVals[idv], na2=80, phiedge=true)
Dcomp[idv]=norm(Dval-DataGS,Inf)
end
return Dcomp
end


function makeqmctab(Dcomp,NVals,NxVals)
printf(fmt::String, args...) = @eval @printf($fmt, $(args...))
levels=length(NxVals)
Nrange=length(NVals)
printf("%s :","Nx\\N")
for np=1:Nrange
    printf("%10d   ", NVals[np])
end
    printf("\n")
for lp=1:levels
    printf("%4d :  ", NxVals[lp])
    for np=1:Nrange
    printf("%8.5e  ", Dcomp[lp,np])
    end
    printf("\n")
end
end


function validate(N, levels, DataGS, s=1.0)
# Get the QMC results
Nx = 50
Diffs=zeros(levels,)
for il=0:levels-1
DataQMC=ctk_qmc_test(N, Nx; s=s, plotme=false, tol=1.e-8)
EDiff=(DataGS-DataQMC)./DataGS;
Diffs[il+1]=norm(EDiff,Inf);
#EDiff=(DataGS-DataQMC)
#Diffs[il+1]=norm(EDiff,1)/norm(DataGS,1)
Nx *= 2
end
return Diffs
end

function readdata(s)
DataGS=zeros(11,2);
s == 1.0 ? sstr=string(1) : sstr="Inf"
fname=string("Siewert:s=",sstr,".dat")
readtab(fname,DataGS)
return DataGS
end
