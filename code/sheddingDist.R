library(fitdistrplus); library(mixtools)

dat <- read.csv('file:///C:/Users/jcanevari/Documents/PhD/Writings/PhD_papers/4.Cb_Shedding/code/yieldOnCbStatus/data/daily_prod_dat(V03).csv', header = TRUE, strip.white = TRUE)

tdat <- data.frame(cowKey = unique(dat$CowKey))
tdat$cc.cb <- dat$cc.cb[match(tdat$cowKey, dat$CowKey)]
tdat$com1.status <- dat$com1.status[match(tdat$cowKey, dat$CowKey)]
tdat$farm <- dat$Farm[match(tdat$cowKey, dat$CowKey)]
tdat$season <- dat$Season[match(tdat$cowKey, dat$CowKey)]
tdat$seroTot <- dat$sero.tot[match(tdat$cowKey, dat$CowKey)]
tdat$parity <- dat$grouped_lact[match(tdat$cowKey, dat$CowKey)]
dat <- tdat

hist(log10(dat$cc.cb[dat$cc.cb>0]), freq = F)

k.con <- log10(dat$cc.cb[dat$cc.cb>0])
k.con <- k.con[!is.na(k.con)]
hist(k.con)

fg <- fitdist(k.con, "gamma")
fln <- fitdist(k.con, "lnorm")
kfw <- fitdist(k.con, "weibull")

plot.legend <- c("gamma","lognormal", "weibull")
denscomp(list(fg, fln, kfw), legendtext = plot.legend)

summary(fln)     

rlnorm(100, meanlog = fln$estimate[1], sdlog = fln$estimate[2])

h.shed <-  k.con[k.con>=3]
l.shed <-  k.con[k.con<3]

median(h.shed)
median(l.shed)
epi.descriptives(l.shed)

#do ratio of median of h.shed with median of l.shed and quantile 25/75 and use that as limits for your prior of epsilon_fL

#alternatively, use the observed distribution and see what happens.

temp <- 10^rlnorm(100, meanlog = fln$estimate[1], sdlog = fln$estimate[2])/10^quantile(k.con, .90)
length(temp[temp>1])


#compare density distributions between shedding at kidding and abortion

dat <- read.csv('C:\\Temp\\abortions_pcr(V02).csv', header = TRUE)

head(dat)
dat$Result[dat$Result == 'rep'] <- 'pos'
dat <- dat[dat$Result == 'pos',]

library(epiR)
epi.descriptives(dat$Calc.conc..copies.uL.)

a.con <- log10(dat$Calc.conc..copies.uL.)
epi.descriptives(a.con)
# n     mean       sd      q25      q50      q75   lower    upper      min      max na
# 14 5.824594 2.283954 3.556715 6.129367 7.623228 4.51539 7.133798 2.245513 8.421604  0
hist(a.con,10)

#-----------------------------------------
# try putting together abortions and on time kiddings and fit mixed distribution
temp <- c(a.con, k.con)
hist(temp)
mixmdl <- normalmixEM(temp)
plot(mixmdl,which=2)


#Do your own plot.

hist(temp, main='',xlab = expression(paste(italic("C. burnetii"), " GE /",  mu, "L", " (log10)")), freq = F, ylim = c(0,.5), col='gray')
lines(seq(0,10,length.out = 100),dnorm(seq(0,10,length.out = 100), mean = 6.753013, sd = 1.4133656)*0.2482477, lwd=2, lty = 5)
lines(seq(0,10,length.out = 100),dnorm(seq(0,10,length.out = 100), mean = 1.973731, sd = 0.6792576)*0.7517523, lwd=2)
legend('topright', legend=c('low-shedders', 'high-shedders'), lty=c(1,5))



# $lambda
# [1] 0.2482477 0.7517523
# $mu
# [1] 6.753013 1.973731
# $sigma
# [1] 1.4133656 0.6792576

lines(seq(0,10,length.out = 100),dnorm(seq(0,10,length.out = 100), mean = 6.753013, sd = 1.4133656)/4)
lines(seq(0,10,length.out = 100),dnorm(seq(0,10,length.out = 100), mean = 1.973731, sd = 0.6792576)*3/4)

rval <- vector('numeric', 1000)
for(i in 1:1000) {
omega <- 0.75
temp <- runif(1,0,1)
if(temp <= omega){
  rval[i] <- kshed <- rnorm(10000, mean = 1.973731, sd = 0.6792576) 
} else {
  rval [i]<- rnorm(10000,  mean = 6.753013, sd = 1.4133656)
}
}
hist(rval, add=T, freq = F)


#what should be the value of omega
pos <- 14
neg <- (73 - 14)
dat <- as.matrix(cbind(pos, neg))
round(epi.conf(dat, ctype = "prop.single"), digits = 3)



#-----------------------------------------

library(ks);library(logspline)
lsfit <- logspline(k.con, lbound=0,ubound =10)
rval <- rlogspline(100000,lsfit)
hist(rval, freq=F)
kfw <- fitdist(rval, "weibull") #fitting for kiddings of infected
kfw <- fitdist(rval, "lnorm") #fitting for kiddings of infected


lsfit <- logspline(a.con, lbound=0,ubound =10)
rval <- rlogspline(100000,lsfit)
epi.descriptives(rval)
hist(rval, freq=F)
afw <- fitdist(rval, "weibull") #fitting for abortions

plot(density(k.con))
lines(density(a.con))

summary(kfw)
summary(afw)
temp <- seq(0,10,length.out =100)
lines(temp,dweibull(temp, shape = 2.869776, scale = 6.454916 ), lty=2) #abortions
lines(temp,dweibull(temp, shape = 1.582719, scale = 2.939030 ), lty=2) #kiddings

qweibull(.90, shape = 2.869776, scale = 6.454916) #abortions; this will become 1
qweibull(.95, shape = 1.582719, scale = 2.939030) #kiddings


#====================================================
# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))

# Draw the boxplot and the histogram 
par(mar=c(0, 4, 1.1, 2.1))
boxplot(a.con , horizontal=TRUE , ylim=c(0,10), xaxt="n" , frame=F)
par(mar=c(4, 4, 1.1, 2.1))
hist(log10(dat$cc.cb[dat$cc.cb>0]), freq = F, main='', col = 'gray', border = F, xlab = expression(paste(italic("C. burnetii"), " GE /",  mu, "L", " (log10)")))
lines(temp,dweibull(temp, shape = 2.869776, scale = 6.454916 ), lty=2) #abortions
# lines(temp,dweibull(temp, shape = 1.582719, scale = 2.939030 ), lty=1) #kiddings
lines(temp,dlnorm(temp, meanlog = 0.7805830, sdlog = 0.5726633), lty=1) #kiddings
