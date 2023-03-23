function averop(p,n)
d=2.0*ones(n)/4.0;
d[1]=2.0/3.0;
d[n]=2.0/3.0;
du=ones(n-1)/4.0;
du[1]=1.0/3.0;
dl=ones(n-1)/4.0;
dl[n-1]=1.0/3.0;
T=Tridiagonal(dl,d,du);
A=sparse(T^p);
return A
end
