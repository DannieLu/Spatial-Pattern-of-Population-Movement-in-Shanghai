######################################## Loading Libraries ######################################
#library(maps)
library(akima)
library(geoR)
library(xtable)
library(MASS)

######################################### Read Data #############################################
opar <- par()
data <- read.csv('Spatial_Project/GSM_all_withDist.csv',header = T)
data <- data[, -1]

################################### Exploratory Data Analysis ###################################
int.data <- interp(data$long, data$lat, data$freq)
par(mfrow=c(1,1))
image(int.data,xlab = 'longitude',ylab = 'latitude',xlim = range(D$long), ylim = range(D$lat))
contour(int.data,add=TRUE)
persp(int.data,xlab = 'longitude',ylab = 'latitude', zlab = 'temperature',theta = 60, phi = 30)

gdata<-as.geodata(D,coords.col=1:2,data.col=3)

# descriptive stats
summary(D$long)
summary(D$lat)
summary(D$freq)
longitude<-obj[,1];latitude<-obj[,2];temp<-obj[,3]
par(mfrow = c(2,2))
hist(temp,xlab='temperature',main='Histogram of Temperature')
boxplot(temp,xlab='temperature',main='Boxplot of Temperature')
stem(temp)
qqnorm(temp,xlab='Normal Quantiles')
qqline(temp, lwd = 2)
qqplot(runif(300,250,270),temp,xlab='Unifrom Quantiles',
       ylab='Sample Quantiles',main='Unifrom Q-Q Plot')
qqline(temp, distribution = function(p) qunif(p, 250,270),
       prob = c(0.1, 0.6), lwd = 2)
summary(temp)
diff(range(temp))

################################### Check spatial Trend ######################################### 
plot(gdata)

#points(gdata,xlab="Coord X",ylab="Coord Y")
points(gdata,xlab="Coord X",ylab="Coord Y",pt.divide="rank.prop",main='rank')
#points(gdata,xlab="Coord X",ylab="Coord Y",cex.max=1.7,col=gray(seq(1,0.1,l=100)),pt.divide="equal")
points(gdata,xlab="Coord X",ylab="Coord Y",pt.divide="quartile",main='quartile')
points(gdata,xlab="Coord X",ylab="Coord Y",pt.divide="quintile",main='quintile')
points(gdata,xlab="Coord X",ylab="Coord Y",pt.divide="decile",main='decile')

# row & col boxplots
par(mfrow=c(1,2))
boxplot(matrix(obj[,3],nrow = length(unique(obj[,2])), ncol = length(unique(obj[,1]))),
        xlab = 'x', main = 'Boxplot across x coord')
boxplot(matrix(obj[,3],nrow = length(unique(obj[,1])), ncol = length(unique(obj[,2])),
               byrow = T), xlab = 'y', main = 'Boxplot across y coord')

# regression
fit1 <- lm(gdata$data ~ gdata$coords[,1] + gdata$coords[,2])
summary(fit1)
fit2 <- lm(gdata$data ~ gdata$coords[,1] + I(gdata$coords[,1]^2) + gdata$coords[,2]
           + I(gdata$coords[,2]^2))
summary(fit2)
fit3 <- lm(gdata$data ~ gdata$coords[,1] + I(gdata$coords[,1]^2) + gdata$coords[,2]
           + I(gdata$coords[,2]^2) + gdata$coords[,1]*gdata$coords[,2])
summary(fit3)
gres <- gdata
gres$data <- residuals(fit3)
plot(gres)

# boxcox transformation
par(mfrow=c(1,2))
boxcox(fit2,lambda = seq(-10,0,0.1))
boxcox(fit3,lambda = seq(-10,0,0.1))

ldata <- log(gdata$data)
lfit1 <- lm(ldata ~ gdata$coords[,1] + gdata$coords[,2])
summary(lfit1)
lfit2 <- lm(ldata ~ gdata$coords[,1] + I(gdata$coords[,1]^2) + gdata$coords[,2]
            + I(gdata$coords[,2]^2))
summary(lfit2)
lfit3 <- lm(ldata ~ gdata$coords[,1] + I(gdata$coords[,1]^2) + gdata$coords[,2]
            + I(gdata$coords[,2]^2) + gdata$coords[,1]*gdata$coords[,2])
summary(lfit3)
lgres <- gdata
lgres$data <- residuals(lfit3)
plot(lgres)

sum <- rbind(c(summary(fit1)$adj.r.square, summary(fit2)$adj.r.square,
               summary(fit3)$adj.r.square),c(summary(lfit1)$adj.r.square,
                                             summary(lfit2)$adj.r.square, summary(lfit3)$adj.r.square))
sum <- data.frame(sum)
rownames(sum) <- c('Data','Log(Data)')
colnames(sum) <- paste('Model',1:3,sep='')
xtable(sum,digits = 4)
xtable(summary(fit3))

# semivariogram
par(mfrow=c(1,2))
vario.class <- variog(gres,estimator.type="classical")
vario.modul <- variog(gres,estimator.type="modulus")
plot(vario.class,main='Variogram for Residual')
points(vario.modul$u,vario.modul$v,col="red")
legend('bottomright',legend = c('classic','Cressie'),,pch = 1,col = c(1,2))
log.vario.class <- variog(lgres,estimator.type="classical")
log.vario.modul <- variog(lgres,estimator.type="modulus")
plot(log.vario.class,main='Variogram for log(residual)')
points(log.vario.modul$u,log.vario.modul$v,col="red")
legend('bottomright',legend = c('classic','Cressie'),,pch = 1,col = c(1,2))
# correlogram
c1 <- correlog(x = gres$coords[,1], y = gres$coords[,2], z = gres$data, increment = 5, 
               resamp = 100, quiet = TRUE)
plot(c1, main = 'Correlogram for Residuals')
c2 <- correlog(x = gres$coords[,1], y = gres$coords[,2], z = lgres$data, increment = 5, 
               resamp = 100, quiet = TRUE)
plot(c2, main = 'Correlogram for Log(Residuals)')

#################################### Check Anistropy ###########################################
# ESC
source('esc.r')
par(mfrow=c(1,1))
ESC(gres$coords[,1],gres$coords[,2],gres$data)
# directional semivariogram
par(opar)
dir.vario <- variog4(gres)
plot(dir.vario)
# change coordinates
coords <- gres$coords
aniso.coords <- coords.aniso(coords,aniso.pars=c(0,3),reverse=FALSE)
obj.aniso <- cbind(gres$data,aniso.coords)
gres.aniso <- as.geodata(obj.aniso,coords.col=2:3,data.col=1)
# check after transformation
gres.var <- variog4(gres.aniso)
par(mfrow=c(1,2))
plot(gres.var)
ESC(gres.aniso$coords[,1],gres.aniso$coords[,2],gres.aniso$data)

#################################### Model Fitting #############################################
# empirical semivariogram
par(mfrow=c(1,1))
ani.vario.class <- variog(gres.aniso,estimator.type="classical")
ani.vario.modul <- variog(gres.aniso,estimator.type="modulus")
plot(ani.vario.class,main='Variogram with Transformed Coords')
points(ani.vario.modul$u,ani.vario.modul$v,col="red")
legend('bottomright',legend = c('classic','Cressie'),,pch = 1,col = c(1,2))

vario <- ani.vario.class
# Spherical covariance
variofit1 <- variofit(vario,ini.cov.pars=c(1,1000),cov.model="spherical")
variofit1.cressie <- variofit(vario,ini.cov.pars=c(1,1000),cov.model="spherical",
                              weights="cressie")
s1 <- c(variofit1$nugget, variofit1$cov.pars[1], variofit1$cov.pars[2], 
        variofit1$practicalRange, variofit1$value)
s2 <- c(variofit1.cressie$nugget, variofit1.cressie$cov.pars[1], variofit1.cressie$cov.pars[2], 
        variofit1.cressie$practicalRange, variofit1.cressie$value)

# exponential
variofit2 <- variofit(vario,ini.cov.pars=c(1,1000),cov.model="exponential")
variofit2.cressie <- variofit(vario,ini.cov.pars=c(1,1000),cov.model="exponential",
                              weights="cressie")

e1 <- c(variofit2$nugget, variofit2$cov.pars[1], variofit2$cov.pars[2], 
        variofit2$practicalRange, variofit2$value)
e2 <- c(variofit2.cressie$nugget, variofit2.cressie$cov.pars[1], variofit2.cressie$cov.pars[2], 
        variofit2.cressie$practicalRange, variofit2.cressie$value)

# Gaussian covariance
variofit3 <- variofit(vario,ini.cov.pars=c(1,1000),cov.model="gaussian")
variofit3.cressie <- variofit(vario,ini.cov.pars=c(1,1000),cov.model="gaussian",
                              weights="cressie")
g1 <- c(variofit3$nugget, variofit3$cov.pars[1], variofit3$cov.pars[2], 
        variofit3$practicalRange, variofit3$value)
g2 <- c(variofit3.cressie$nugget, variofit3.cressie$cov.pars[1], variofit3.cressie$cov.pars[2], 
        variofit3.cressie$practicalRange, variofit3.cressie$value)

# matern covariance
variofit4 <- variofit(vario,ini.cov.pars=c(1,1000),cov.model="matern",
                      fix.kappa=FALSE,kappa=5)
variofit4.cressie <- variofit(vario,ini.cov.pars=c(1,1000),cov.model="matern",
                              weights="cressie",fix.kappa=FALSE,kappa=5)

m1 <- c(variofit4$nugget, variofit4$cov.pars[1], variofit4$cov.pars[2], 
        variofit4$practicalRange, variofit4$value)
m2 <- c(variofit4.cressie$nugget, variofit4.cressie$cov.pars[1], variofit4.cressie$cov.pars[2], 
        variofit4.cressie$practicalRange, variofit4.cressie$value)

# comparison
com <- rbind(s1,s2,e1,e2,g1,g2,m1,m2)
com <- data.frame(com)
xtable(com, digits = 4)

# Likelihood Functions
# spherical
spherical.ml.fit <- likfit(gres.aniso,cov.model="spherical",ini.cov.pars=c(1,1000))
print(summary(spherical.ml.fit))
spherical.anis.ml.fit <- likfit(gres,cov.model="spherical",ini.cov.pars=c(1,1000),
                                fix.psiA=FALSE,psiA=0,fix.psiR=FALSE,psiR=5)
print(summary(spherical.anis.ml.fit))
s <- c(summary(spherical.ml.fit)$likelihood[3], summary(spherical.ml.fit)$likelihood[4],
       summary(spherical.anis.ml.fit)$likelihood[3], summary(spherical.anis.ml.fit)$likelihood[4])
s <- as.numeric(s)

# exponential
exponential.ml.fit <- likfit(gres.aniso,cov.model="exponential",ini.cov.pars=c(1,1000))
print(summary(exponential.ml.fit))
exponential.anis.ml.fit <- likfit(gres,cov.model="exponential",ini.cov.pars=c(1,1000),
                                  fix.psiA=FALSE,psiA=0,fix.psiR=FALSE,psiR=3)
print(summary(exponential.anis.ml.fit))
e <- c(summary(exponential.ml.fit)$likelihood[3], summary(exponential.ml.fit)$likelihood[4],
       summary(exponential.anis.ml.fit)$likelihood[3], summary(exponential.anis.ml.fit)$likelihood[4])
e <- as.numeric(e)

# Gaussian
gaussian.ml.fit <- likfit(gres.aniso,cov.model="gaussian",ini.cov.pars=c(1,50))
print(summary(gaussian.ml.fit))
gaussian.anis.ml.fit <- likfit(gres,cov.model="gaussian",ini.cov.pars=c(1,50),
                               fix.psiA=FALSE,psiA=0,fix.psiR=FALSE,psiR=5)
print(summary(gaussian.anis.ml.fit))
g <- c(summary(gaussian.ml.fit)$likelihood[3], summary(gaussian.ml.fit)$likelihood[4],
       summary(gaussian.anis.ml.fit)$likelihood[3], summary(gaussian.anis.ml.fit)$likelihood[4])
g <- as.numeric(g)

# Matern
matern.ml.fit <- likfit(gres.aniso,cov.model="matern",ini.cov.pars=c(1,1000))
print(summary(matern.ml.fit))
matern.anis.ml.fit <- likfit(gres,cov.model="matern",ini.cov.pars=c(1,1000),
                             fix.psiA=FALSE,psiA=0,fix.psiR=FALSE,psiR=5)
print(summary(matern.anis.ml.fit))
m <- c(summary(matern.ml.fit)$likelihood[3], summary(matern.ml.fit)$likelihood[4],
       summary(matern.anis.ml.fit)$likelihood[3], summary(matern.anis.ml.fit)$likelihood[4])
m <- as.numeric(m)

# comparison
comp <- rbind(s, e, g, m)
comp <- data.frame(comp)
rownames(comp) <- c('Spherical', 'Exponential', 'Gaussian', 'Matern')
colnames(comp) <- c('isoAIC', 'isoBIC', 'anisoAIC', 'anisoBIC')
xtable(comp, digits = 4)

plot(vario)
lines(variofit1)
lines(variofit1.cressie,col="blue")
lines(gaussian.ml.fit,col="red")
lines(spherical.anis.ml.fit,col="green3")
legend('bottomright',legend=c('Sph.WLS','Sph.cressie','Gau.ML','Sph.ML'),lty=1,
       col=c('black','blue','red','green3'))
#################################### Model Prediction ###########################################
int.gres.aniso <- interp(gres.aniso$coords[,1], gres.aniso$coords[,2], gres.aniso$data)
image(int.gres.aniso)
contour(int.gres.aniso, add = T)
loci.ani <- expand.grid(seq(min(gres.aniso$coords[,1]),max(gres.aniso$coords[,1]),10),
                        seq(min(gres.aniso$coords[,2]),max(gres.aniso$coords[,2]),10))
new <- c(258,44)
newlat <- new[2]*pi/180
newlong <- new[1]*pi/180
new.x <- dx*(newlong-mean.long)/(max.long-min.long)
new.y <- dy*(newlat-mean.lat)/(max.lat-min.lat)
newx <- c(1,new.x, new.x^2,new.y,new.y^2,new.x*new.y)
aniso.coords <- coords.aniso(t(c(new.x,new.y)),aniso.pars=c(0,3),reverse=FALSE)
# Gaussian Prediction 
# Image
gauss.pred.image <- krige.conv(gres.aniso,locations=loci.ani,
                               krige=krige.control(type.krige="OK",obj.model=gaussian.ml.fit))

image(gauss.pred.image)
contour(gauss.pred.image,add=TRUE)

# Value
gauss.pred <- krige.conv(gres.aniso,locations=aniso.coords,
                         krige=krige.control(type.krige="OK",obj.model=gaussian.ml.fit))
gauss.pred$predict
gauss.pred$krige.var
gauss.pred$predict + sum(fit3$coefficients*newx)
gauss.pred$predict + sum(fit3$coefficients*newx)-1.96*sqrt(gauss.pred$krige.var)
gauss.pred$predict + sum(fit3$coefficients*newx)+1.96*sqrt(gauss.pred$krige.var)

# Spherical Semivariogram Prediction 
# Image
sph.pred.image <- krige.conv(gres.aniso,locations=loci.ani,
                             krige=krige.control(type.krige="OK",obj.model=variofit1))

image(sph.pred.image)
contour(sph.pred.image,add=TRUE)

# Value
sph.pred <- krige.conv(gres.aniso,locations=aniso.coords,
                       krige=krige.control(type.krige="OK",obj.model=variofit1))
sph.pred$predict
sph.pred$krige.var
sph.pred$predict + sum(fit3$coefficients*newx)
sph.pred$predict + sum(fit3$coefficients*newx)-1.96*sqrt(sph.pred$krige.var)
sph.pred$predict + sum(fit3$coefficients*newx)+1.96*sqrt(sph.pred$krige.var)
