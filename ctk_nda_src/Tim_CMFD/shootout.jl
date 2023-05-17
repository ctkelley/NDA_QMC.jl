function shootout()
Ndims=[4, 8, 16, 32, 32]
Nxdims=[512, 512, 2048, 1024, 2048]
Ncdims=[64, 64, 128, 128, 128]
nfd=length(Ndims)
for i=1:nfd
figure(i)
plot_test(4096*Ndims[i], Nxdims[i], Ncdims[i])
end
end
