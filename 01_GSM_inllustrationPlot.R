# GSM spatial
# semivariogram plot for each hour

setwd('C:/Users/Dannie/Dropbox/Dissertation/GMSData/R/SpatialAnalysis')
setwd('C:/Users/dlu/Dropbox/Dissertation/GMSData/R/SpatialAnalysis')

library(RandomFields) # for vtti R version 
library(geoR)
library(akima)
library(gstat)
library(ncf) # for correlation 

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
semivar.plot2<-function(D,log.transf,xmax)
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
    simul.var <- variog(simul.geo,estimator.type="classical",breaks=seq(0,xmax,0.0005))
    #simul.var <- variog(simul.geo,estimator.type="modulus",breaks=seq(0,0.4,0.01))
    Var.tab<-rbind(Var.tab,simul.var$v)
  }
  # plot variance
  limx<-c(0,xmax)
  limy<-ceiling(max(Var.tab)*1.05)
  for (i in 1:length(uniq.hour))
  {
    plot(simul.var$u,Var.tab[i,],pch=18,cex=0.7,col=i,xlim=limx,ylim=c(0,limy),yaxt='n',xaxt='n',xlab=' ',ylab=' ')
    par(new=T)
    print(i)
  }
  axis(side=2,at=seq(0,limy,length.out=5),tcl=0.4,lwd.ticks=1,mgp=c(0,0.2,0))
  print(limy)
  axis(side=1,at=seq(0,xmax,length.out=5),tcl=0.4,lwd.ticks=1,mgp=c(0,0.2,0))
  mtext(side=1,text="Coordination Distance",line = 1.5)
  mtext(side=2,text="Semivariance",line=1.5)
  legend('topleft',bty='n',pch=rep(18,length(uniq.hour)), col=seq(1,length(uniq.hour),1),
         legend = c('5:00-6:00','6:00-7:00','7:00-8:00','8:00-9:00','9:00-10:00'),cex=0.8)
  par(new=F)
}

# 4. LSE fit lines
LSE.fit<-function(D,Model,log.transf,xmax)
# D: dataset;  
# D$freq.k in 1000 scale
# Model: covariance model
# log.transf==1: log transformation £» ==0: original scale
# ini.cov : initial partial sill(sigma^2) & range (phi)
# output: 1 Plot with fitted lines for each hour
{
  D$lgfreq<-log(D$freq+1)
  D$freq.k<-D$freq/1000
  uniq.hour<-unique(D$hour)
  limx<-c(0,xmax)
  paratab<-data.frame()
  for (i in uniq.hour)
  {
    #simul.geo <- as.geodata(D[D$hour==i,c('long','lat','lgfreq')])
    if (log.transf==1)
    {simul.geo <- as.geodata(D[D$hour==i,c('long','lat','lgfreq')]); limy=25}
    if (log.transf==0)
    {simul.geo <- as.geodata(D[D$hour==i,c('long','lat','freq.k')]);limy=300 }
    simul.var <- variog(simul.geo,estimator.type="classical",breaks=seq(0,xmax,0.0005),max.dist=xmax)
    #simul.modvar <- variog(simul.geo,estimator.type="modulus",breaks=seq(0,0.4,0.01))
    plot(simul.var,pch=18,cex=0.7,col=i-4,xlim=limx,ylim=c(0,limy),yaxt='n',xaxt='n',xlab='',ylab='')

    sill<-max(simul.var$v)*0.9
    fit.gau <- variofit(simul.var,ini.cov.pars=c(sill,0.1),cov.model=Model)
    lines(fit.gau,col=i-4,lwd=2)  
    par(new=T)
    print(i)
    print(fit.gau)
    para<-round(c(fit.gau$nugget,fit.gau$cov.pars,fit.gau$practicalRange,fit.gau$value),4)
    paratab<-rbind(paratab,para)
  }
  axis(side=2,at=seq(0,limy,length.out=5),tcl=0.4,lwd.ticks=1,mgp=c(0,0.2,0))
  axis(side=1,at=seq(0,xmax,length.out=5),tcl=0.4,lwd.ticks=1,mgp=c(0,0.2,0))
  mtext(side=1,text="Coordination Distance",line = 1.5)
  mtext(side=2,text="Semivariance",line=1.5)
  mtext(side=3,text=Model,line=0.5)
  legend('topleft',bty='n',pch=rep(18,length(uniq.hour)),lty=rep('solid',length(uniq.hour)), col=1:5,
         legend = c('5:00-6:00','6:00-7:00','7:00-8:00','8:00-9:00','9:00-10:00'),cex=1)
  par(new=F)
  
  names(paratab)<-c('nugget','sigma.sq','phi','range','sum.of.sq')
  return(paratab)
}




Corr<-function(D,StepLen,xmax)
  # plot correlation
  # D: dataset;  
  # D$freq.k in 1000 scale
  # output: 6th row: mean of each bin
{
  D$freq.k<-D$freq/1000
  uniq.hour<-unique(D$hour)
  limx<-c(0,xmax)
  corrtab<-data.frame()
  for (i in uniq.hour)
  {
    Dsub<-D[D$hour==i,]
    print(i)
    ncf.cor <- correlog(Dsub$long, Dsub$lat,Dsub$freq.k, increment=0.01,resamp = 100)
    plot(ncf.cor$mean.of.class,ncf.cor$correlation,type='l',
         pch=18,cex=0.7,col=i-4,xlim=limx,ylim=c(-0.3,1),yaxt='n',xaxt='n',xlab='',ylab='')
    par(new=T)
    corrtab<-rbind(corrtab,ncf.cor$correlation)
  }
  corrtab<-rbind(corrtab,ncf.cor$mean.of.class)
  axis(side=2,at=round(seq(-0.2,1,length.out=7),1),tcl=0.4,lwd.ticks=1,mgp=c(0,0.2,0))
  axis(side=1,at=round(seq(0,xmax,length.out=5),2),tcl=0.4,lwd.ticks=1,mgp=c(0,0.2,0))
  mtext(side=1,text="Coordination Distance",line = 1.5)
  mtext(side=2,text="Correlation",line=1.5)
  mtext(side=3,text=" ",line=0.5)
  legend(x=0.085,y=0.9,bty='n',pch=rep(18,5),lty=rep('solid',5), col=1:5,
         legend = c('5:00-6:00','6:00-7:00','7:00-8:00','8:00-9:00','9:00-10:00'),cex=1)
  par(new=F) 
  return(corrtab)
    
}


Cortab<-Corr(D,StepLen=0.01,xmax=0.15)
write.csv(Cortab,'Corrtab.csv')

Cortab<-read.csv('Corrtab.csv')

xmax=0.15
ind<-Cortab[6,]<=xmax

limx=c(0,xmax)

for( i in 1:5)
{
  x<-as.numeric(Cortab[6,ind])
  y<-as.numeric(Cortab[i,ind])
  plot(x,y,type = 'l',pch=18,cex=0.7,col=i,xlim=limx,ylim=c(-0.3,1),yaxt='n',xaxt='n',xlab=' ',ylab=' ')
  par(new=T)
}
axis(side=2,at=round(seq(-0.2,1,length.out=7),1),tcl=0.4,lwd.ticks=1,mgp=c(0,0.2,0))
axis(side=1,at=round(seq(0,xmax,length.out=5),2),tcl=0.4,lwd.ticks=1,mgp=c(0,0.2,0))
mtext(side=1,text="Coordination Distance",line = 1.5)
mtext(side=2,text="Correlation",line=1.5)
mtext(side=3,text=" ",line=0.5)
legend(x=0.09,y=1,bty='n',lty=rep('solid',5), col=1:5,
       legend = c('5:00-6:00','6:00-7:00','7:00-8:00','8:00-9:00','9:00-10:00'),cex=0.9,y.intersp = 0.7)
par(new=F) 


##################
# Plot
hist.plot(D)

semivar.plot(D)

semivar.plot2(D,1,0.4) # log
semivar.plot2(D,0,xmax = 0.15)

#log scale
LSE.fit(D,'gaussian',1,c(10,0.3))  # log
LSE.fit(D,'matern',1,c(20,0.3))
LSE.fit(D,'exponential',1,c(10,0.3))

#orginal 
LSE.fit(D,'gaussian',0,xmax=0.4)  # orginal 
Gau<-LSE.fit(D,'gaussian',0,xmax=0.15)  # orginal
Cub<-LSE.fit(D,'cubic',0,xmax=0.15) 
Sph<-LSE.fit(D,'spherical',0,xmax=0.15) 
Mat<-LSE.fit(D,'matern',0,xmax=0.15) 
Exp<-LSE.fit(D,'exponential',0,xmax=0.15) 

#############
D$lgfreq<-log(D$freq+1)

  simul.geo <- as.geodata(D[,c('long','lat','lgfreq')])
  #simul.geo <- as.geodata(D[,c('long','lat','freq')])
  
  simul.var <- variog(simul.geo,estimator.type="classical",breaks=seq(0,0.4,0.0005))
  #simul.modvar <- variog(simul.geo,estimator.type="modulus",breaks=seq(0,0.4,0.0005))
  plot(simul.var,pch=18,cex=0.8) 
  #plot(simul.modvar,pch=18,cex=0.8)


#############

