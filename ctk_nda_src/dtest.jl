function dtest(n)
dx=1.0/(n-1); x=.5*dx:dx:1.0-.5*dx; se=sin.(4.0*pi*x); 
D1s=spdiagm(n-2, n-1, 0 => -ones(n-2,), 1 => ones(n-2,))
D1s ./= dx;
ca=D1s*se
xi=dx:dx:1-dx; ce=4.0*pi*cos.(4.0*pi*xi);
println(norm(ce-ca,Inf))
end
