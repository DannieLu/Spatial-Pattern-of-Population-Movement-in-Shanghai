# GSM spatial
# semivariogram plot for each hour

setwd('C:/Users/Dannie/Dropbox/Dissertation/GMSData/R/SpatialAnalysis')

library(geoR)
library(akima)
library(gstat)


###################
# start from here
D<-read.csv('GSM_all.csv')

########################
# plot functions

# 1. histogram of freq and log freq for each hour  5:00 ~ 9:00 
hist.plot<-function(D)
{
  uniq.hour<-unique(D$hour)
  par(mfcol=c(2,length(uniq.hour)),mar=c(4,4,2,1))
    for (i in uniq.hour)
  {
     hist(D[D$hour==i,'freq'],main=' ',xlab="freq")
    hist( log(D[D$hour==i,'freq']), main='',xlab='logfreq')
  }
  par(mfcol=c(1,1))
}

# 2. semivarigram plot: 
#    10 sub plots: 5 hour, 2 method: classical and cressie
semivar.plot<-function(D)
{
  D$lgfreq<-log(D$freq+1)
  uniq.hour<-unique(D$hour)
  par(mfcol=c(2,length(uniq.hour)),mar=c(4,4,2,1))
  for (i in uniq.hour)
  {
    #simul.geo <- as.geodata(D[D$hour==i,c('long','lat','lgfreq')])
    simul.geo <- as.geodata(D[D$hour==i,c('long','lat','freq')])
    simul.var <- variog(simul.geo,estimator.type="classical",breaks=seq(0,0.4,0.01))
    simul.modvar <- variog(simul.geo,estimator.type="modulus",breaks=seq(0,0.4,0.01))
    plot(simul.var,pch=18,cex=0.8) 
    plot(simul.modvar,pch=18,cex=0.8)
  }
  par(mfcol=c(1,1))
}

# 3. semivarigram plot: 
#    1 plot: different colors for each hour; method: classical 
semivar.plot2<-function(D,log.transf)
  # log.transf==1: log transformation
  #          ==0: original scale
{
  D$lgfreq<-log(D$freq+1)
  uniq.hour<-unique(D$hour)
  Var.tab<-data.frame()
  for (i in uniq.hour)
  {
    if (log.transf==1)
    {simul.geo <- as.geodata(D[D$hour==i,c('long','lat','lgfreq')])}
    if (log.transf==0)
    {simul.geo <- as.geodata(D[D$hour==i,c('long','lat','freq')])  }
    simul.var <- variog(simul.geo,estimator.type="classical",breaks=seq(0,0.4,0.01))
    #simul.var <- variog(simul.geo,estimator.type="modulus",breaks=seq(0,0.4,0.01))
    Var.tab<-rbind(Var.tab,simul.var$v)
  }
  # plot variance
  limx<-c(0,0.4)
  limy<-ceiling(max(Var.tab)*1.05)
  for (i in 1:length(uniq.hour))
  {
    plot(simul.var$u,Var.tab[i,],pch=18,cex=0.7,col=i,xlim=limx,ylim=c(0,limy),yaxt='n',xaxt='n',xlab=' ',ylab=' ')
    par(new=T)
    print(i)
  }
  axis(side=2,at=seq(0,limy,length.out=5),tcl=0.4,lwd.ticks=1,mgp=c(0,0.2,0))
  print(limy)
  axis(side=1,at=seq(0,0.4,0.1),tcl=0.4,lwd.ticks=1,mgp=c(0,0.2,0))
  mtext(side=1,text="Coordination Distance",line = 1.5)
  mtext(side=2,text="Semivariance",line=1.5)
  legend('topleft',bty='n',pch=rep(18,length(uniq.hour)), col=seq(1,length(uniq.hour),1),
         legend = c('5:00-6:00','6:00-7:00','7:00-8:00','8:00-9:00','9:00-10:00'),cex=0.7)
  par(new=F)
}

# 4. LSE fit lines
LSE.fit<-function(D,Model,log.transf,ini.cov)
# D: dataset;  Model: covariance model
# ini.cov : partial sill(sigma^2) & range (phi)
# output: 1 Plot with fitted lines for each hour
{
  D$lgfreq<-log(D$freq+1)
  uniq.hour<-unique(D$hour)
  limx<-c(0,0.4)
  for (i in uniq.hour)
  {
    #simul.geo <- as.geodata(D[D$hour==i,c('long','lat','lgfreq')])
    if (log.transf==1)
    {simul.geo <- as.geodata(D[D$hour==i,c('long','lat','lgfreq')]); limy=25}
    if (log.transf==0)
    {simul.geo <- as.geodata(D[D$hour==i,c('long','lat','freq')]);limy=3*10^8 }
    simul.var <- variog(simul.geo,estimator.type="classical",breaks=seq(0,0.4,0.01))
    #simul.modvar <- variog(simul.geo,estimator.type="modulus",breaks=seq(0,0.4,0.01))
    fit.gau <- variofit(simul.var,ini.cov.pars=ini.cov,cov.model=Model)
    plot(simul.var,pch=18,cex=0.7,col=i,xlim=limx,ylim=c(0,limy),yaxt='n',xaxt='n',xlab='',ylab='')
    lines(fit.gau,col=i)  
    par(new=T)
    print(i)
  }
  axis(side=2,at=seq(0,limy,5),tcl=0.4,lwd.ticks=1,mgp=c(0,0.2,0))
  axis(side=1,at=seq(0,0.4,0.1),tcl=0.4,lwd.ticks=1,mgp=c(0,0.2,0))
  mtext(side=1,text="Coordination Distance",line = 1.5)
  mtext(side=2,text="Semivariance",line=1.5)
  par(new=F)
  
}

##################
# Plot
hist.plot(D)

semivar.plot(D)

semivar.plot2(D,1) # log
semivar.plot2(D,0)

LSE.fit(D,'gaussian',1,c(10,0.3))  
LSE.fit(D,'gaussian',0,c(10,0.3))  
LSE.fit(D,'matern',1,c(20,0.3))
LSE.fit(D,'exponential',1,c(10,0.3))


