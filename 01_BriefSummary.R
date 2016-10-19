# GSM Data Spatial Analysis
# 

setwd('C:/Users/Dannie/Documents/GitHub/SpatialProject')


library(geoR)
library(akima)
library(gstat)


D.raw<-read.table('fishnetpoint09.txt',sep=',',header = T)
D<-D.raw[c('long','lat','freq')]

summary(D$long)
summary(D$lat)
summary(D$freq)

plot(D$long,D$lat)
boxplot(D$freq)
hist(D$freq,breaks=20,main='',xlab = 'Number of mobile user',ylab = 'Frequency of sites')

lgfreq<-log(D$freq+1)
hist(lgfreq)
D$lgfreq<-lgfreq

# contour and persp
int.scp <- interp(D$long, D$lat, D$lgfreq)
image(int.scp,xlim=range(D$long),ylim=range(D$lat))

contour(int.scp,add=TRUE)

persp(int.scp,theta=60, phi=30)

# semivariogram
simul.geo <- as.geodata(D[,c(1,2,4)])
simul.var <- variog(simul.geo,estimator.type="classical",breaks=seq(0,0.4,0.01))
simul.modvar <- variog(simul.geo,estimator.type="modulus",breaks=seq(0,0.4,0.01))# (Hawkins and Cressie,

plot(simul.var)
plot(simul.modvar)


# check anisotropy
dir.variog<-variog4(simul.geo, direction = c(0, pi/4, pi/2, 3*pi/4))
plot(dir.variog,xlim=c(0,4),bty='n',cex=0.6,bty='n')


# Variogram fit: LSE
fit.nug.cressie.exp <- variofit(simul.var,ini.cov.pars=c(10,0.3),cov.model="exponential",weights="cressie",fix.kappa = F)
fit.nug.cressie.mat <- variofit(simul.var,ini.cov.pars=c(10,0.3),cov.model='matern',weights="cressie",fix.kappa = F)
fit.nug.cressie.gau <- variofit(simul.var,ini.cov.pars=c(10,0.3),cov.model="gaussian",weights="cressie",fix.kappa = F)
fit.nug.cressie.sph <- variofit(simul.var,ini.cov.pars=c(10,0.3),cov.model="spherical",weights="cressie",fix.kappa = F)
fit2 <- variofit(simul.var,ini.cov.pars=c(10,0.3),cov.model="cubic",weights="cressie",fix.kappa = F)

fit.nug.cressie.exp <- variofit(simul.modvar,ini.cov.pars=c(10,0.3),cov.model="exponential",weights="cressie",fix.kappa = F)
fit.nug.cressie.mat <- variofit(simul.modvar,ini.cov.pars=c(10,0.3),cov.model='matern',weights="cressie",fix.kappa = F)
fit.nug.cressie.gau <- variofit(simul.modvar,ini.cov.pars=c(10,0.3),cov.model="gaussian",weights="cressie",fix.kappa = F)
fit.nug.cressie.sph <- variofit(simul.modvar,ini.cov.pars=c(10,0.3),cov.model="spherical",weights="cressie",fix.kappa = F)
fit2 <- variofit(simul.modvar,ini.cov.pars=c(10,0.3),cov.model="cubic",weights="cressie",fix.kappa = F)


plot(simul.var,main='',pch=20)
lines(fit.nug.cressie.exp,col='black')
lines(fit.nug.cressie.mat,col='red')
lines(fit.nug.cressie.gau,col='blue')
lines(fit.nug.cressie.sph,col='green')
lines(fit2,col='orange')
legend('topleft',col = c('black','red','blue','green','orange'),bty='n',lty=1,cex = 0.75,
       legend = c('exponential','matern','gaussian','spherical','cubic'))

summary(fit2)

#other model "circular", "cubic", "wave","power", "powered.exponential"
fit1 <- variofit(simul.var,ini.cov.pars=c(12,0.3),cov.model="circular",weights="cressie",fix.kappa = F)
fit2 <- variofit(simul.var,ini.cov.pars=c(12,0.3),cov.model="cubic",weights="cressie",fix.kappa = F)
fit3<- variofit(simul.var,ini.cov.pars=c(12,0.3),cov.model="power",weights="cressie",fix.kappa = F)
fit4 <- variofit(simul.var,ini.cov.pars=c(12,0.3),cov.model="powered.exponential",weights="cressie",fix.kappa = F)

plot(simul.var,main='')
lines(fit1,col='black')
lines(fit2,col='red')
lines(fit3,col='blue')
lines(fit4,col='green')
legend('topleft',col = c('black','red','blue','green'),bty='n',lty=1,cex = 0.75,
       legend = c('exponential','matern','gaussian','spherical'))


# Max likelihood method
#--------------------------------------------------------
ini=c(10,0.3)
ml.exp <- likfit(simul.geo,cov.model="exponential",ini.cov.pars=ini,fix.kappa = F)
ml.mat <- likfit(simul.geo,cov.model="matern",ini.cov.pars=ini,fix.kappa = F)
ml.gau <- likfit(simul.geo,cov.model="gaussian",ini.cov.pars=c(5,0.2))
ml.sph <- likfit(simul.geo,cov.model="spherical",ini.cov.pars=ini,fix.kappa = F)
ml.cub <- likfit(simul.geo,cov.model="cubic",ini.cov.pars=ini,fix.kappa = F)

plot(simul.var,main='',pch=20)
lines(ml.exp,col='black')
lines(ml.mat,col='red')
lines(ml.gau,col='blue')
lines(ml.sph,col='green')
lines(ml.cub,col='orange')
legend('topleft',col = c('black','red','blue','green','orange'),bty='n',lty=1,cex = 0.75,
       legend = c('exponential','matern','gaussian','spherical','cubic'))


summary(fit.nug.cressie.gau )
