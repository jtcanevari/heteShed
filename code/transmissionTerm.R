SNP=1
beta=2.16
E=1
SNP*(1-beta*2.718^-E)


E=seq(0,4,1)
propInf <- SNP*(1-(2.718)^-E)

plot(E, y= propInf)
