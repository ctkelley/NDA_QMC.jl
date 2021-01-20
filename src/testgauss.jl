function testgauss(m)
x, w = gausslegendre(m)
# Map to [0,1]
w .*= .5; x .+= 1.0; x .*= .5;
errs=zeros(2*m,)
for p=0:2*m-1
    exact=1.0/(1.0+p)
    y=x.^p
#    approx=y'*w
    approx=dot(w,y)
    errs[p+1]=abs(exact-approx)
end
maximum(errs)
end
