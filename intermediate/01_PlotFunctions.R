
######################################## Loading Libraries ######################################
#library(maps)
library(akima)
library(geoR)
library(xtable)
library(MASS)
library(gstat)
library(ncf) # for correlation 
library(dplyr)

# This file writes functions for many types of plots of variable frequency for each hour 5:00 ~ 9:00

########################
# exploratory plot functions

# 1. Image Plots
image.combine <- function(D){
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    int.data <- interp(Di$long, Di$lat, Di$freq)
    image(int.data, xlab = paste(i,': 00 ~',i+1,': 00 am'), ylab = '',
          xlim = range(Di$long), ylim = range(Di$lat), zlim = range(D$freq))
    contour(int.data, add=TRUE)
  }
  par(mfcol=c(1, 1))
}

# 2. Persp Plots
persp.combine <- function(D){
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    int.data <- interp(Di$long, Di$lat, Di$freq)
    persp(int.data, xlab = paste(i,': 00 ~',i+1,': 00 am'), ylab = '', 
          zlab = 'population density (pop/sq.km)',
          xlim = range(Di$long), ylim = range(Di$lat), zlim = range(D$freq))
  }
  par(mfcol=c(1,1))
}
  
# 3. histograms  
hist.combine<-function(D)
{
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    fr <- as.vector(select(Di, freq)$freq)
    if (max(fr) > 10000){
      br <- seq(min(fr), max(fr) + 1000, 1000)
    } else {
      br<- seq(0,15,0.5)
    }
    hist(fr, breaks=br, col=alpha("black", 0.5), border = 'white',
        # xlim = limx, ylim = c(0,limy),
        xlab = 'population density (pop/sq.km)',
        ylab = 'frequency',main = paste(i,': 00 ~',i+1,': 00 am'))
  }
  par(mfcol=c(1,1))
}

# 4. Boxplots
boxplot.combine <- function(D){
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    fr <- as.vector(select(Di, freq)$freq)
    boxplot(fr, main = paste(i,': 00 ~',i+1,': 00 am'))
  }
  par(mfcol=c(1,1))
}

# 5. QQplots
qq.combine <- function(D){
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    fr <- as.vector(select(Di, freq)$freq)
    qqnorm(fr, main = paste(i,': 00 ~',i+1,': 00 am'))
    qqline(fr, lwd = 2)
  }
  par(mfcol=c(1,1))
}

####################
# spatial trend plot functions

# 6. row & col boxplots
long.combine <- function(D){
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    boxplot(freq ~ long, data = Di, xlab = 'longitude', main = paste(i,': 00 ~',i+1,': 00 am'),
            names = round(as.vector(arrange(distinct(Di, long), long)$long), digits = 2))
  }
  par(mfcol=c(1,1))
}
lat.combine <- function(D){
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    boxplot(freq ~ lat, data = Di, xlab = 'latitude', main = paste(i,': 00 ~',i+1,': 00 am'),
            names = round(as.vector(arrange(distinct(Di, lat), lat)$lat), digits = 2))
  }
  par(mfcol=c(1,1))
}

########################
# semivariograms

# 7. empirical semivarigram plots
#    5 sub plots: 5 hour, 2 method: classical and cressie
semivar.combine<-function(D)
{
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    gDi <- as.geodata(Di, coords.col = 1:2, data.col = 3)
    vario.class <- variog(gDi, estimator.type="classical", breaks=seq(0,0.4,0.01))
    vario.modul <- variog(gDi, estimator.type="modulus", breaks=seq(0,0.4,0.01))
    plot(vario.class, main = paste(i,': 00 ~',i+1,': 00 am'))
    points(vario.modul$u, vario.modul$v, col="red")
    legend('bottomright', legend = c('classic','Cressie'),pch = 1, cex = 0.8, col = c(1,2))
  }
  par(mfcol=c(1,1))
}

#    1 plots: different colors for each hour; method: Cressie weight 
semivar.inone<-function(D, xmax)
{
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  Var.tab <- data.frame()
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    gDi <- as.geodata(Di, coords.col = 1:2, data.col = 3)
    vario.modul <- variog(gDi, estimator.type="modulus", breaks=seq(0,0.4,0.01))
    Var.tab <- rbind(Var.tab, vario.modul$v)
  }
  # plot variance
  limx<-c(0,xmax)
  limy<-ceiling(max(Var.tab)*1.05)
  for (i in 1:length(uniq.hour))
  {
    plot(vario.modul$u, Var.tab[i,], pch = 18, cex = 0.7, col = i, xlim = limx, ylim = c(0,limy),
         yaxt = 'n', xaxt = 'n', xlab = ' ', ylab = ' ')
    par(new = T)
    print(i)
  }
  axis(side = 2, at = seq(0, limy, length.out = 5), tcl = 0.4, lwd.ticks = 1, mgp = c(0,0.2,0))
  print(limy)
  axis(side = 1, at = seq(0, xmax, length.out = 5), tcl = 0.4, lwd.ticks = 1, mgp = c(0,0.2,0))
  mtext(side = 1, text = "Coordination Distance", line = 1.5)
  mtext(side = 2, text = "Semivariance", line = 1.5)
  legend('topleft', bty = 'n', pch = 18, col = 1:length(uniq.hour),
         legend = c('5:00-6:00','6:00-7:00','7:00-8:00','8:00-9:00','9:00-10:00'), cex=0.8)
  par(new = F)
}

###################
# check anisotropy

# 8. ESC plot
esc.combine <- function(D){
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    source('esc.r')
    ESC(Di$long,Di$lat, Di$fre)
  }
  par(mfcol=c(1,1))
}

# 9. directional semivariogram
dir.combine <- function(D){
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    gDi <- as.geodata(Di, coords.col = 1:2, data.col = 3)
    dir.vario <- variog4(gDi)
    plot(dir.vario, xlab = paste(i,': 00 ~',i+1,': 00 am'))
  }
  par(mfcol=c(1,1))
}

# 4. LSE fit lines
LSE.fit<-function(D,Model,log.transf,xmax)
# D: dataset;  
# D$freq.k in 1000 scale
# Model: covariance model
# log.transf==1: log transformation ?? ==0: original scale
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


