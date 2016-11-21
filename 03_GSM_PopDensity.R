# GSM Data Process
# population density distribution over time


setwd('C:/Users/Dannie/Dropbox/Dissertation/GMSData/R/SpatialAnalysis')
setwd('C:/Users/dlu/Dropbox/Dissertation/GMSData/R/SpatialAnalysis')

D<-read.csv('GSM_all.csv')

hist(D$freq)
boxplot(D$freq)
summary(D$freq)

##########
# hist
Hist.plot(D,Log='F',limx = c(0,80000) )
Hist.plot(D,Log='F',limx = c(0,60000) )
Hist.plot(D,Log='T',limx = c(0,12) )

# combined plot
Gamma.plot(D,limx=80000,plothist = 'T') #with all histgrams
Gamma.plot(D,limx=20000,plothist = 'T') #with all histgrams
Gamma.plot(D,limx=80000,plothist = 'F') # with overall histogram
Gamma.plot(D,limx=20000,plothist = 'F') # with overall histogram

# seperate plots
Gamma.plot2(D,limx=20000)

# Log scale
Logfit(D,dist='gamma',limx0=3,limx=13)
Logfit(D,dist='weibull',limx0=3,limx=13)

#-----FUNCTIONS-----------------------------------------------------------------
library(scales) # alpha function for transparent color

Hist.plot<-function(D,limx,Log)
{
    trun<-0
    par(mfrow=c(2,3),mar=c(4,4,2,1))
    #limy<-0.0002
    limy<-160
    for (i in 5:9)
    {
      if (Log=='T')
        { den<- log(D[D$hour==i & D$freq>=trun,'freq']+1); br<- seq(0,20,0.5) }
      if (Log=='F') 
        { den<-D[D$hour==i & D$freq>=censor,'freq']; br<- seq(min(den),max(den)+1000,1000) }
      
      hist(den,xlim = limx,ylim = c(0,limy),breaks=br,
           col=alpha("black",0.5),border = 'white',xlab = " population density (pop/sq.km)",ylab = 'frequency',main = '')
      legend('topright',bty='n',cex=1,legend = paste(i,':00~',i+1,':00 am'))
    }
    par(mfcol=c(1,1))
}

####Plot Gamma and hist for each hour in one plot
Gamma.plot<-function(D,limx,plothist)
{
  cols=c('red2','gold','forestgreen','blue','black')
  if(plothist=='F')
  {
    hist(D[D$freq!=0,"freq"],xlim = c(0,limx),breaks=seq(min(D$freq),max(D$freq)+1000,1000),col=alpha('black',0.3),border = 'white',
         freq = F,xlab = "Population Density",main='Gamma Distribution')
    limy<-0.00015
    for (i in 5:9)
    {
      den<-D[D$hour==i & D$freq>0,'freq']
      shape = mean(den)^2/var(den)  
      scale = var(den)/mean(den) 
      x<-seq(0,limx,0.01)
      par(new=T)
      curve(dgamma(x,shape=shape,scale = scale),from = 0,to=limx,ylim=c(0,limy),col=cols[i-4],xlab = '',ylab='',xaxt='n',yaxt='n')
      print(paste('hour=',i));print(paste('shape=',shape));print(paste('scale=',scale))
      if (i<9){ par(new=T)}
    }
    
  }# # end: if(plothist=='F')
 
 if(plothist=='T')
 {
   limy<-0.0002
   for (i in 5:9)
   {
     den<-D[D$hour==i & D$freq>0,'freq']
     hist(den,xlim = c(0,limx),ylim = c(0,limy),breaks=seq(min(den),max(den)+1000,1000),freq = F,
          col=alpha(i,0.3),border = 'white',xlab = " ",ylab = '')
     shape = mean(den)^2/var(den)  
     scale = var(den)/mean(den) 
     x<-seq(0,limx,0.01)
     par(new=T)
     curve(dgamma(x,shape=shape,scale = scale),from = 0,to=limx,ylim=c(0,limy),col=cols[i-4],xlab = '',ylab='',xaxt='n',yaxt='n')
     print(paste('hour=',i));print(paste('shape=',shape));print(paste('scale=',scale))
     if (i<9){ par(new=T)}
   }
 } # end: if(plothist=='T')
  legend('topright',bty='n',lty='solid', col=cols,
         legend = c('5:00-6:00','6:00-7:00','7:00-8:00','8:00-9:00','9:00-10:00'),cex=0.7) 
}
# end: gamma.plot function

####Plot hist and Gamma for each hour 
Gamma.plot2<-function(D,limx)
{
  for (i in 5:9)
  {
    den<-D[D$hour==i,'freq']
    hist(den,xlim = c(0,limx),breaks=seq(min(den),max(den)+1000,1000),col=alpha('black',0.3),border = 'white',
         freq = F,xlab = "Population Density",main='Gamma Distribution')
    shape = mean(den)^2/var(den)  
    scale = var(den)/mean(den) 
    x<-seq(0,limx,0.01)
    par(new=T)
    curve(dgamma(x,shape=shape,scale = scale),from = 0,to=limx,col='red',xlab = '',ylab='',xaxt='n',yaxt='n')
    legend('topright',legend = paste(i,':00~',i+1,':00'),bty='n',cex=0.8)
    print(paste('hour=',i));print(paste('shape=',shape));print(paste('scale=',scale))
    
  }
}
# end :Gamma.plot2



####Plot log density: Gamma for each hour in one plot
Logfit<-function(D,dist,limx0,limx)
  # dist: gamma; weibull
{
  cols=c('red2','gold','forestgreen','blue','black')
  den<-log(D[ D$freq!=0,'freq'])
  hist(den,xlim = c(limx0,limx),breaks=seq(min(den),max(den)+0.5,0.5),
       xlab = "Log Population Density",main=dist,col=alpha('black',0.2),border = 'white')
  for (i in 5:9)
  {
    den<-log(D[D$hour==i & D$freq!=0,'freq']+1)
    if (dist=='gamma')
    {
      shape = mean(den)^2/var(den)  
      scale = var(den)/mean(den) 
      x<-seq(limx0,limx,0.01)
      par(new=T)
      curve(dgamma(x,shape=shape,scale = scale),from = limx0,to=limx,col=cols[i-4],xlab = '',ylab='',xaxt='n',yaxt='n')
    }
    if (dist=='weibull')
    {
      pa1<-fitdistr(den,'weibull')
      shape = as.numeric(pa1$estimate[1]  )
      scale = as.numeric(pa1$estimate[2])
      x<-seq(limx0,limx,0.01)
      par(new=T)
      curve(dweibull(x,shape=shape,scale = scale),from =limx0,to=limx,col=cols[i-4],xlab = '',ylab='',xaxt='n',yaxt='n') 
    }
    print(paste('hour=',i));print(paste('shape=',shape));print(paste('scale=',scale))
    if (i<9){ par(new=T)}
  }
  legend('topright',bty='n',lty='solid', col=cols,
         legend = c('5:00-6:00','6:00-7:00','7:00-8:00','8:00-9:00','9:00-10:00'),cex=0.7)
  
}





###########################################################################################################
#-----------------------------------
limx0<-5
limx<-15
den<-log(D[ D$freq!=0,'freq'])
hist(den,xlim = c(limx0,limx),breaks=seq(min(den),max(den)+0.5,0.5),
     xlab = "Population Density",main='Gamma Distribution',col=alpha('black',0.3),border = 'white')
for (i in 5:9)
{
  den<-log(D[D$hour==i & D$freq!=0,'freq']+1)
  shape = mean(den)^2/var(den)  
  scale = var(den)/mean(den) 
  x<-seq(0,limx,0.01)
  par(new=T)
 curve(dgamma(x,shape=shape,scale = scale),from = limx0,to=limx,col=i,xlab = '',ylab='',xaxt='n',yaxt='n')
  #curve(dweibull(x,shape=shape,scale = scale),from = 0,to=limx,col=i,xlab = '',ylab='',xaxt='n',yaxt='n')
  print(paste('hour=',i));print(paste('shape=',shape));print(paste('scale=',scale))
  if (i<9){ par(new=T)}
}
legend('topright',bty='n',lty='solid', col=seq(5,9,1),
       legend = c('5:00-6:00','6:00-7:00','7:00-8:00','8:00-9:00','9:00-10:00'),cex=0.7)



#--------------------------------------
library(MASS)

wei<-log(D[ D$freq!=0,'freq'])

pa1<-fitdistr(wei,'weibull')
hist(wei,xlim=c(0,14))
par(new=T)
curve(dweibull(x,shape=as.numeric(pa1$estimate[1]),scale = as.numeric(pa1$estimate[2])),from = 0,to=14,xlab = '',ylab='',xaxt='n',yaxt='n')



fit_weibull <- function(x)
{
  xbar <- mean(x)
  varx <- var(x)
  f <- function(b){return(gamma(1+2/b)/gamma(1+1/b)^2 - 1 - varx/xbar^2)}
  bhat <- uniroot(f,c(0.02,50))$root
  ahat <- xbar/gamma(1+1/bhat)
  return(c(ahat,bhat))
}

pa.wei<-fit_weibull(log(D[ D$freq!=0,'freq']))
hist(log(D[ D$freq!=0,'freq']),xlim=c(0,14))
par(new=T)
curve(dweibull(x,shape=pa.wei[1],scale = pa.wei[2]),from = 0,to=14,xlab = '',ylab='',xaxt='n',yaxt='n')
