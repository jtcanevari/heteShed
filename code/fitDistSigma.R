#fit a gamma distribution to sigma and use that in the model

setwd('C:\\Users\\jcanevari\\Documents\\Projects\\PhD\\qFeverABM\\data\\abc\\monoExp4')

files = list.files(pattern="*.csv")
myfiles = lapply(files, read.csv)

parNames <- c("gamma","muE", "epsilon_p","phi", "fI", "iP","alpha","rho","seed")

for(i in 1:length(myfiles)){
  colnames(myfiles[[i]])=parNames
}

df <- do.call("rbind", myfiles)
head(df)

hist(df$iP)

library(fitdistrplus )

fit.gamma <- fitdist(df$iP, distr = "gamma", method = "mle")
summary(fit.gamma)
plot(fit.gamma)

scale = 2.911

hist(rgamma(1000, shape = 15.34, scale=2.911))
