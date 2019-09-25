library(Rcpp); library(BH)
Sys.getenv("PKG_CXXFLAGS")
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
#------------------------------------------------------------------------------
#source models
#------------------------------------------------------------------------------

setwd("~/Documents/Projects/PhD/heteShed/code")
sourceCpp("qF_mixDistList.cpp")
sourceCpp("qF_mixDistArray.cpp")
#------------------------------------------------------------------------------
#choose parameter values
#------------------------------------------------------------------------------
par <- c('gamma' = 0.05 , 
         'mu' = 0.2 ,
         'phi' = 0.85 ,
         'fI' = 0.3,
         'iP' =  45,
         'alpha' = 0.3 ,
         'omega' = 0.2  ,
         'seed' = runif(1) )

#------------------------------------------------------------------------------
#run model
#------------------------------------------------------------------------------
library(dplyr)

#run list version
out <- qFmixDistList(par) 
plot(out$time, out$IP, type='l')

#run array version
out <- as.data.frame(qFmixDistArray(par))
colnames(out) <- c('year','parity', 'nkid','nkif', 'nabo') 

rval<- out %>% group_by(year) %>% summarise_each(list(sum)) #grouping by year but not parity
rval$kidPrev <- rval$nkif/rval$nkid
plot(rval$year, rval$kidPrev, type = 'o')

rval$aboPrev <- rval$nabo/(rval$nkid+rval$nabo)
lines(rval$year, rval$aboPrev, type = 'l',col = 'red')



