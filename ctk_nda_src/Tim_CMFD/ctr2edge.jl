function ctr2edge(cval)
N=length(cval)
eval=zeros(N+1)
eval[2:N]=.5*(cval[1:N-1]+cval[2:N])
eval[1]=2.0*cval[1]-eval[2]
eval[N+1]=2.0*cval[N]-eval[N]
#tst1=1.5*cval[1]-.5*cval[2]
#tstN=1.5*cval[N]-.5*cval[N-1]
#println(tst1-eval[1],"  ",tstN-eval[N+1])
return eval
end

function xftest(N)
h=1.0/N; z=collect(.5*h:h:1.0-.5*h);
x=collect(0.0:h:1.0);
sc=sin.(5.0*pi*z);
se=ctr2edge(sc);
sol=sin.(5.0*pi*x);
println(norm(sol-se,Inf))
end
