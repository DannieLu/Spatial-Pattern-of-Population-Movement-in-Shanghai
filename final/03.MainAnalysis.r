######################################## Loading Libraries ######################################
#library(maps)
library(akima)
library(geoR)
library(xtable)
library(MASS)
#library(ggplot2)
library(dplyr)
#library(gstat)
#library(ncf) # for correlation 
library(scales)
library(mgcv)

##################################### Data Manipulation #########################################
# 00 file

######################################### Read Data #############################################
opar <- par()
data <- read.csv('Spatial_Project/GSM_all_withDist.csv',header = T)
data <- data[, -1]
#summary(data)
xtable(summary(data))

######################################### Functions #############################################
source('Spatial_Project/intermediate/01_PlotFunctions.R')

################################### Exploratory Data Analysis ###################################
image.combine(data)
persp.combine(data)

# descriptive stats
hist.combine(data)
lgdata <- data
lgdata$freq <- log(data$freq + 1)
hist.combine(lgdata)
boxplot.combine(data)
boxplot.combine(lgdata)
qq.combine(data)
qq.combine(lgdata)

################################### Check spatial Trend ######################################### 
gdata5<-as.geodata(filter(data, hour == 5), coords.col=2:3,data.col=4)
plot(gdata5)
gdata6<-as.geodata(filter(data, hour == 6), coords.col=2:3,data.col=4)
plot(gdata6)
gdata7<-as.geodata(filter(data, hour == 7), coords.col=2:3,data.col=4)
plot(gdata7)
gdata8<-as.geodata(filter(data, hour == 8), coords.col=2:3,data.col=4)
plot(gdata8)
gdata9<-as.geodata(filter(data, hour == 9), coords.col=2:3,data.col=4)
plot(gdata9)

#points(gdata,xlab="Coord X",ylab="Coord Y")
#points(gdata,xlab="Coord X",ylab="Coord Y",pt.divide="rank.prop",main='rank')
#points(gdata,xlab="Coord X",ylab="Coord Y",cex.max=1.7,col=gray(seq(1,0.1,l=100)),pt.divide="equal")
#points(gdata,xlab="Coord X",ylab="Coord Y",pt.divide="quartile",main='quartile')
#points(gdata,xlab="Coord X",ylab="Coord Y",pt.divide="quintile",main='quintile')
#points(gdata,xlab="Coord X",ylab="Coord Y",pt.divide="decile",main='decile')

# row & col boxplots
long.combine(data)
lat.combine(data)

# regression
fit1 <- lm(data ~ .^2 + I(long^2) + I(lat^2), data = gdata5)
xtable(summary(fit1))
stepAIC(fit1)
fit2 <- lm(data ~ .^2 + I(long^2) + I(lat^2), data = gdata6)
xtable(summary(fit2))
stepAIC(fit2)
fit3 <- lm(data ~ .^2 + I(long^2) + I(lat^2), data = gdata7)
xtable(summary(fit3))
stepAIC(fit3)
fit4 <- lm(data ~ .^2 + I(long^2) + I(lat^2), data = gdata8)
xtable(summary(fit4))
stepAIC(fit4)
fit5 <- lm(data ~ .^2 + I(long^2) + I(lat^2), data = gdata9)
xtable(summary(fit5))
stepAIC(fit5)

# boxcox transformation
par(mfrow=c(2,3))
boxcox(lm(data + 1 ~ .^2 + I(long^2) + I(lat^2), data = gdata5))
boxcox(lm(data + 1 ~ .^2 + I(long^2) + I(lat^2), data = gdata6))
boxcox(lm(data + 1 ~ .^2 + I(long^2) + I(lat^2), data = gdata7))
boxcox(lm(data + 1 ~ .^2 + I(long^2) + I(lat^2), data = gdata8))
boxcox(lm(data + 1 ~ .^2 + I(long^2) + I(lat^2), data = gdata9))
par(mfrow=c(1,1))

lgdata5<-as.geodata(filter(lgdata, hour == 5), coords.col=2:3,data.col=4)
plot(lgdata5)
lgdata6<-as.geodata(filter(lgdata, hour == 6), coords.col=2:3,data.col=4)
plot(lgdata6)
lgdata7<-as.geodata(filter(lgdata, hour == 7), coords.col=2:3,data.col=4)
plot(lgdata7)
lgdata8<-as.geodata(filter(lgdata, hour == 8), coords.col=2:3,data.col=4)
plot(lgdata8)
lgdata9<-as.geodata(filter(lgdata, hour == 9), coords.col=2:3,data.col=4)
plot(lgdata9)
# row & col boxplots
long.combine(lgdata)
lat.combine(lgdata)

# regression
lfit1 <- lm(data ~ .^2 + I(long^2) + I(lat^2), data = gdata5)
summary(lfit1)
stepAIC(lfit1)
lfit2 <- lm(data ~ .^2 + I(long^2) + I(lat^2), data = gdata6)
summary(lfit2)
stepAIC(lfit2)
lfit3 <- lm(data ~ .^2 + I(long^2) + I(lat^2), data = gdata7)
summary(lfit3)
stepAIC(lfit3)
lfit4 <- lm(data ~ .^2 + I(long^2) + I(lat^2), data = gdata8)
summary(lfit4)
stepAIC(lfit4)
lfit5 <- lm(data ~ .^2 + I(long^2) + I(lat^2), data = gdata9)
summary(lfit5)
stepAIC(lfit5)

gres5 <- gdata5
gres6 <- gdata6
gres7 <- gdata7
gres8 <- gdata8
gres9 <- gdata9
gres5$data <- fit1$residuals
gres6$data <- fit2$residuals
gres7$data <- fit3$residuals
gres8$data <- fit4$residuals
gres9$data <- fit5$residuals

plot(gres5)
plot(gres6)
plot(gres7)
plot(gres8)
plot(gres9)

gres5 <- data.frame(cbind(gres5$coords, gres5$data))
gres5$hour <- 5
names(gres5)[3] <- 'fre'
gres6 <- data.frame(cbind(gres6$coords, gres6$data))
gres6$hour <- 6
names(gres6)[3] <- 'fre'
gres7 <- data.frame(cbind(gres7$coords, gres7$data))
gres7$hour <- 7
names(gres7)[3] <- 'fre'
gres8 <- data.frame(cbind(gres8$coords, gres8$data))
gres8$hour <- 8
names(gres8)[3] <- 'fre'
gres9 <- data.frame(cbind(gres9$coords, gres9$data))
gres9$hour <- 9
names(gres9)[3] <- 'fre'
gres <- bind_rows(gres5, gres6, gres7, gres8, gres9)

#lgres5 <- gdata5
#lgres6 <- gdata6
#lgres7 <- gdata7
#lgres8 <- gdata8
#lgres9 <- gdata9
#lgres5$data <- lfit1$residuals
#lgres6$data <- lfit2$residuals
#lgres7$data <- lfit3$residuals
#lgres8$data <- lfit4$residuals
#lgres9$data <- lfit5$residuals

#plot(lgres5)
#plot(lgres6)
#plot(lgres7)
#plot(lgres8)
#plot(lgres9)

#################################### Check Anistropy ###########################################
# ESC
esc.combine(gres)

# directional semivariogram
dir.combine(gres)

#################################### Model Fitting #############################################
# semivariogram
semivar.combine(gres)
semivar.inone(gres, 0.4)
semivar.inone(gres, 0.4/sqrt(2))

# empirical semivariogram fitting
semfit_s <- semi.fit(gres, Model = 'spherical', xmax = 0.15)
xtable(semfit_s)
semfit_e <- semi.fit(gres, Model = 'exponential', xmax = 0.15)
xtable(semfit_e)
semfit_g <- semi.fit(gres, Model = 'gaussian', xmax = 0.15)
xtable(semfit_g)
semfit_c <- semi.fit(gres, Model = 'cubic', xmax = 0.15)
xtable(semfit_c)
semfit_m <- semi.fit(gres, Model = 'matern', xmax = 0.15)
xtable(semfit_m)
semfit_ci <- semi.fit(gres, Model = 'circular', xmax = 0.15)
xtable(semfit_ci)
semfit_p <- semi.fit(gres, Model = 'power', xmax = 0.15)
xtable(semfit_p)
semfit_pe <- semi.fit(gres, Model = 'powered.exponential', xmax = 0.15)
xtable(semfit_pe)

# comparison
com <- rbind(semfit_s[, 5], semfit_e[, 5], semfit_g[, 5], semfit_c[, 5], semfit_m[, 5],
             semfit_ci[, 5], semfit_p[, 5], semfit_pe[, 5])
me <- apply(com, 1, mean)
com <- cbind(com, me)
rm(me)
com <- data.frame(com)
colnames(com) <- c(5:9, 'mean')
rownames(com) <- c('Spherical', 'Exponential', 'Gaussian', 'Cubic', 'Matern', 'Circular',
                   'Power', 'Powered.exponential')
which.min(com$mean)
xtable(com, digits = 4)
semfit_s

# Likelihood fitting
likefit_s <- lik.fit(gres, Model = 'spherical', xmax = 0.15, r = 0.1)
xtable(likefit_s)
likefit_e <- lik.fit(gres, Model = 'exponential', xmax = 0.15, r = 0.1)
xtable(likefit_e)
likefit_g <- lik.fit(gres, Model = 'gaussian', xmax = 0.15, r = 0.001)
xtable(likefit_g)
likefit_c <- lik.fit(gres, Model = 'cubic', xmax = 0.15, r = 0.1)
xtable(likefit_c)
likefit_m <- lik.fit(gres, Model = 'matern', xmax = 0.15, r = 0.1)
xtable(likefit_m)

# comparison
comp <- rbind(likefit_s[, 5], likefit_e[, 5], likefit_g[, 5], likefit_c[, 5], likefit_m[, 5], 
              likefit_s[, 6], likefit_e[, 6], likefit_g[, 6], likefit_c[, 6], likefit_m[, 6])
ma <- apply(comp, 1, mean)
comp <- cbind(comp, ma)
rm(ma)
comp <- data.frame(comp)
colnames(comp) <- c(5:9, 'mean')
rownames(comp) <- c('Spherical AIC', 'Exponential AIC', 'Gaussian AIC', 'Cubic AIC', 'Matern AIC',
                   'Spherical BIC', 'Exponential BIC', 'Gaussian BIC', 'Cubic BIC', 'Matern BIC')
which.min(comp[1:5, 6])
which.min(comp[6:10, 6])
xtable(comp, digits = 4)
likefit_e
likefit_m

#################################### Model Prediction ###########################################
compare <- pred.compare(data, gres)

pred.semi.combine(data, gres)

pred.lik.combine(data, gres)


#################################### GAM with Predictor #########################################
glm_nsp <- glm(freq ~ D2Metro * D2Road + I(D2Metro^2) + I(D2Road^2),
               family = Gamma, data = data[data$freq > 0, ])
xtable(summary(glm_nsp))
stepAIC(glm_nsp, direction = 'both')
sum(glm_nsp$residuals^2)
gam_sp <- gam(freq ~ D2Metro * D2Road + I(D2Metro^2) + I(D2Road^2) + s(long, lat), 
              family = Gamma, data = data[data$freq > 0, ])
xtable(summary(gam_sp)$p.table)
AIC(gam_sp)
gam_sp1 <- gam(freq ~ D2Metro * D2Road + I(D2Metro^2) + s(long, lat), 
              family = Gamma, data = data[data$freq > 0, ])
xtable(summary(gam_sp1)$p.table, digits = 4)
AIC(gam_sp1)


glm_nsp_nb <- glm(freq ~ D2Metro * D2Road + I(D2Metro^2) + I(D2Road^2), family = nb, data = data)
xtable(summary(glm_nsp_nb))
sum(glm_nsp_nb$residuals^2)
gam_sp_nb <- gam(freq ~ D2Metro * D2Road + I(D2Metro^2) + I(D2Road^2) + s(long, lat), family = nb, data = data)
xtable(summary(gam_sp_nb)$p.table)
AIC(gam_sp_nb)


