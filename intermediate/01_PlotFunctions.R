
######################################## Loading Libraries ######################################
#library(maps)
library(akima)
library(geoR)
library(xtable)
library(MASS)
#library(gstat)
#library(ncf) # for correlation 
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
  
# 3. Histograms  
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

#########################
# spatial trend plot functions

# 6. Row & Col boxplots
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

#########################
# semivariograms

# 7. Empirical Semivarigram Plots
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

#########################
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

# 9. Directional Semivariogram
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

#########################
# model fitting

# 10. Empirical Semivariogram Fitting
semi.fit <- function(D, Model, xmax)
{# D: dataset;  
# D$freq.k in 1000 scale
# Model: covariance model
# ini.cov : initial partial sill(sigma^2) & range (phi)
# output: 1 Plot with fitted lines for each hour

  D$freq.k <- D$fre/1000
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  limx <- c(0, xmax)
  paratab <- data.frame()
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    gDi <- as.geodata(select(Di, long, lat, freq.k)); limy = 150
    #simul.var <- variog(gDi, estimator.type = "classical", breaks = seq(0, xmax, 0.0005),
    #                    max.dist = xmax)
    simul.var <- variog(gDi, estimator.type = "modulus", breaks = seq(0, xmax, 0.0005),
                           max.dist = xmax)
    plot(simul.var, pch = 18, cex = 0.7, col = i - 4, xlim = limx, ylim = c(0, limy),
         yaxt = 'n', xaxt = 'n', xlab = '', ylab = '')

    sill <- max(simul.var$v)*0.9
    fit.gau <- variofit(simul.var, ini.cov.pars = c(sill, 0.1), cov.model = Model, 
                        weights = 'cressie')
    lines(fit.gau, col = i-4, lwd = 2)  
    par(new = T)
    #print(i)
    #print(fit.gau)
    para <- round(c(fit.gau$nugget, fit.gau$cov.pars, fit.gau$practicalRange, fit.gau$value), 4)
    paratab <- rbind(paratab, para)
  }
  axis(side = 2, at = seq(0, limy, length.out = 5), tcl = 0.4, lwd.ticks = 1, mgp = c(0, 0.2, 0))
  axis(side = 1, at = seq(0, xmax, length.out = 5), tcl = 0.4, lwd.ticks = 1, mgp = c(0, 0.2, 0))
  mtext(side = 1, text = "Coordination Distance", line = 1.5)
  mtext(side = 2, text = "Semivariance", line = 1.5)
  mtext(side = 3, text = Model, line = 0.5)
  legend('topleft', bty = 'n', pch =18, lty = 'solid', col = 1:5,
         legend = c('5:00-6:00', '6:00-7:00', '7:00-8:00', '8:00-9:00', '9:00-10:00'), cex = 1)
  par(new = F)
  
  names(paratab) <- c('nugget', 'sigma.sq', 'phi', 'range', 'sum.of.sq')
  rownames(paratab) <- 5:9
  return(paratab)
}

# 11. Likelihood Function Fitting
lik.fit <- function(D, Model, xmax, r)
{# D: dataset;  
  # D$freq.k in 1000 scale
  # Model: covariance model
  # ini.cov : initial partial sill(sigma^2) & range (phi)
  # output: 1 Plot with fitted lines for each hour
  
  D$freq.k <- D$fre/1000
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  limx <- c(0, xmax)
  paratab <- data.frame()
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    gDi <- as.geodata(select(Di, long, lat, freq.k)); limy = 150
    #simul.var <- variog(gDi, estimator.type = "classical", breaks = seq(0, xmax, 0.0005),
    #                    max.dist = xmax)
    simul.var <- variog(gDi, estimator.type = "modulus", breaks = seq(0, xmax, 0.0005),
                        max.dist = xmax)
    plot(simul.var, pch = 18, cex = 0.7, col = i - 4, xlim = limx, ylim = c(0, limy),
         yaxt = 'n', xaxt = 'n', xlab = '', ylab = '')
    
    sill <- max(simul.var$v)*0.9
    mlfit <- likfit(gDi, cov.model = Model, ini.cov.pars=c(sill, r))
    lines(mlfit, col = i-4, lwd = 2)  
    par(new = T)
    print(i)
    print(mlfit)
    para <- c(mlfit$nugget, mlfit$cov.pars, mlfit$practicalRange, 
                    summary(mlfit)$likelihood[3:4])
    para <- unlist(para)
    names(para) <- NULL
    para <- round(para, 4)
    paratab <- rbind(paratab, para)
  }
  axis(side = 2, at = seq(0, limy, length.out = 5), tcl = 0.4, lwd.ticks = 1, mgp = c(0, 0.2, 0))
  axis(side = 1, at = seq(0, xmax, length.out = 5), tcl = 0.4, lwd.ticks = 1, mgp = c(0, 0.2, 0))
  mtext(side = 1, text = "Coordination Distance", line = 1.5)
  mtext(side = 2, text = "Semivariance", line = 1.5)
  mtext(side = 3, text = Model, line = 0.5)
  legend('topleft', bty = 'n', pch =18, lty = 'solid', col = 1:5,
         legend = c('5:00-6:00', '6:00-7:00', '7:00-8:00', '8:00-9:00', '9:00-10:00'), cex = 1)
  par(new = F)
  
  names(paratab) <- c('nugget', 'sigma.sq', 'phi', 'range', 'AIC', 'BIC')
  rownames(paratab) <- 5:9
  return(paratab)
}


##########################
# model prediction

# 12. Model Prediction
      # Compare between four models
pred.compare <- function(D)
{
  uniq.hour <- as.vector(distinct(D, hour)$hour)
  op <- par(no.readonly = TRUE)
  par(mfrow = c(2, 2), mar = c(4,3,1,1), mgp = c(1.5, 0.5, 0))
  compare <- data.frame()
  true <- data.frame()
  sph <- data.frame()
  expo <- data.frame()
  mate <- data.frame()
  model <- list(fit1, fit2, fit3, fit4, fit5)
  for (i in uniq.hour)
  {
    Di <- filter(D, hour == i)
    gDi <- as.geodata(Di, coords.col = c('long', 'lat'), data.col = 'freq')
    #loci <- expand.grid(seq(min(gDi$coords[, 1]), max(gDi$coords[, 1]), 0.005),
    #                        seq(min(gDi$coords[, 2]),max(gDi$coords[, 2]), 0.005))
    int.gdata <- interp(gDi$coords[, 1], gDi$coords[, 2], gDi$data)
    loci <- expand.grid(int.gdata$x, int.gdata$y)
    
    # Spherical Semivariogram Prediction Model
    simul.var <- variog(gDi, estimator.type = "modulus")
    sill <- max(simul.var$v)*0.9
    variofit_s <- variofit(simul.var, ini.cov.pars = c(sill, 0.1), cov.model = 'spherical', 
                           weights = 'cressie')
    sph.pred.image <- krige.conv(gDi, locations = loci,
                                 krige = krige.control(type.krige = "OK", obj.model = variofit_s))
    newx <- cbind(1, loci[, 1], loci[, 2], I(loci[, 1]^2), I(loci[, 2]^2), loci[, 1]*loci[, 2])
    pred_s <- sph.pred.image$predict + newx %*% coef(model[[i-4]])
    
    # Exponential Likelihood Prediction Image
    mlfit_e <- likfit(gDi, cov.model = 'exponential', ini.cov.pars=c(sill, 0.1))
    exp.pred.image <- krige.conv(gDi, locations = loci,
                                 krige = krige.control(type.krige = "OK", obj.model = mlfit_e))
    pred_e <- exp.pred.image$predict + newx %*% coef(model[[i-4]])
    
    # Matern Likelihood Prediction Image
    mlfit_m <- likfit(gDi, cov.model = 'matern', ini.cov.pars=c(sill, 0.1))
    matern.pred.image <- krige.conv(gDi, locations = loci,
                                    krige = krige.control(type.krige = "OK", obj.model = mlfit_m))
    pred_m <- matern.pred.image$predict + newx %*% coef(model[[i-4]])
    
    zmax <- max(gDi$data, pred_s, pred_e, pred_m)
    zmin <- min(gDi$data, pred_s, pred_e, pred_m)
    # true image
    plot(loci, type="n", xlab = 'True Image', ylab = '')
    image(int.gdata, add = T, zlim = c(zmin, zmax)) 
    contour(int.gdata, add=TRUE)
    
    # semivariogram spherical
    plot(loci, type="n", xlab = 'Semivariogram Spherical', ylab = '')
    image(sph.pred.image, values = pred_s, add = T, zlim = c(zmin, zmax))
    contour(sph.pred.image, add=TRUE)
    
    # likelihood exponential
    plot(loci, type="n", xlab = 'Likelihood Exponential', ylab = '')
    image(exp.pred.image, values = pred_e, add = T, zlim = c(zmin, zmax))
    contour(exp.pred.image, add=TRUE)
    
    # likelihood matern
    plot(loci, type="n", xlab = 'Likelihood Matern', ylab = '')
    image(matern.pred.image, values = pred_m, add = T, zlim = c(zmin, zmax))
    contour(matern.pred.image, add=TRUE)
    
    par(op)
    
    # compare models using residual sum of squares
    rss_s <- sum((int.gdata$z - matrix(pred_s, 40, 40))^2)
    rss_e <- sum((int.gdata$z - matrix(pred_e, 40, 40))^2)
    rss_m <- sum((int.gdata$z - matrix(pred_m, 40, 40))^2)
    comp <- c(rss_s, rss_e, rss_m)
    compare <- rbind(compare, comp)
    colnames(compare) <- c('Spherical', 'Exponential', 'Matern')
    rownames(compare) <- i
    
    # write out predicted values
    tru <- cbind(rep(int.gdata$x, 40), rep(int.gdata$y, each = 40), as.vector(int.gdata$z), i)
    true <- rbind (true, tru)
    colnames(true) <- c('x', 'y', 'value', 'hour')
    sphe <- cbind(rep(int.gdata$x, 40), rep(int.gdata$y, each = 40), pred_s, i)
    sph <- rbind (sph, sphe)
    expon <- cbind(rep(int.gdata$x, 40), rep(int.gdata$y, each = 40), pred_e, i)
    expo <- rbind (expo, expon)
    mate <- cbind(rep(int.gdata$x, 40), rep(int.gdata$y, each = 40), pred_m, i)
    mater <- rbind (mate, mater)
  }
  write.csv('true', file = 'Spatial_Project/intermediate/true.csv')
  write.csv('sph', file = 'Spatial_Project/intermediate/spherical.csv')
  write.csv('expo', file = 'Spatial_Project/intermediate/exponential.csv')
  write.csv('mate', file = 'Spatial_Project/intermediate/matern.csv')
  return(compare)
}

      # Compare between time
pred.semi.combine <- function(D){
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