
# Areal Data
# GSM Shanghai

# Text book, Page 163, Problem 12

library(RandomFields)
library(geoR)
library(akima)
library(gstat)

library(coda) # for spBayes
library(abind) # for spBayes
library(magic) # for spBayes
library(Formula) # for spBayes
library(spBayes)

library(Matrix)
library(sp)
library(spdep)
library(maps)
library(maptools)
library(classInt)
library(RColorBrewer)


setwd('//vtti.ad.vt.edu/Data/Users/dlu/Private/Documents/GitHub/Spatial_Project/ArealData')

## Boundary & polygon

mn.poly =   readShapePoly("XY1038_1143_to_TAZ_Project.shp")  # map at noon 

# mn.poly =   readShapePoly("XY2359_0457_to_TAZ_Project.shp")  # map at night

plot(mn.poly,asp=1,main='Boundary of Shanghai')

summary(mn.poly$Density)
boxplot(mn.poly$Density)
hist(mn.poly$Density)



## corepleth map

den.break<-seq(0,100000,20000)
den.break1<-c(0,5000,10000,20000,40000,80000,100000)
spplot(mn.poly, "Density", at = den.break1, col.regions = rev(brewer.pal(10,"RdBu")),main='Population Density')


##SAR and CAR modeling using spdep package  

##polygon and neighbour
mn.nb = poly2nb(mn.poly)
summary(mn.nb)
mn.adj.mat = nb2mat(mn.nb, style="B")
table(as.vector(mn.adj.mat))
rn <- sapply(slot(mn.poly, "polygons"), function(x) slot(x, "ID"))


##Compute Moran's I and Geary's C
mn.listw = nb2listw(mn.nb, style="B", zero.policy=TRUE)
print(mn.listw,zero.policy=T)
#mn.listw = nb2listw(mn.nb, style="W", zero.policy=TRUE)

mn.moran.out = moran.test(mn.poly$Density, listw=mn.listw , zero.policy=TRUE)
mn.geary.out = geary.test(mn.poly$Density, listw=mn.listw , zero.policy=TRUE)
print(mn.moran.out)
print(mn.geary.out)


##SAR model for proportion


## without covarites

# errorsarlm  
mn.sar.out = spautolm(Density~ 1, data=mn.poly, family="SAR", listw=mn.listw, zero.policy=TRUE)
summary(mn.sar.out)
mn.sar.fitted = fitted(mn.sar.out)
mn.poly$fitted.sar = mn.sar.fitted

# lagsarlm
mn.lagsar.out = lagsarlm(Density~ 1, data=mn.poly, listw=mn.listw, zero.policy=TRUE)
summary(mn.lagsar.out)



mn.coords = coordinates(mn.poly)
mn.lat<-as.numeric(mn.coords[,1])
mn.lng<-as.numeric(mn.coords[,2])

mn.sar.out = spautolm(data.proportion~ mn.lat+mn.lng, data=mn.poly, family="SAR", listw=mn.listw, zero.policy=TRUE)
summary(mn.sar.out)


##CAR model regressing rates.FT on NWBIR79.FT
mn.car.out = spautolm(Density~ 1, data=mn.poly, family="CAR", listw=mn.listw, zero.policy=TRUE)
summary(mn.car.out)
mn.car.fitted = fitted(mn.car.out)
mn.poly$fitted.car = mn.car.fitted

mn.car.out2 = spautolm(data.proportion~ data.screen, data=mn.poly, family="CAR", listw=mn.listw, zero.policy=TRUE)
summary(mn.car.out2) # higer AIC, no good


##Draw the maps using the maps function

# color
color.pallete = rev(brewer.pal(6,"RdBu"))

class.raw = classIntervals(var=mn.poly$Density, n=6, style="fixed", fixedBreaks=den.break1, dataPrecision=5)
color.code.raw = findColours(class.raw, color.pallete)

class.fitted.car = classIntervals(var=mn.poly$fitted.car, n=6, style="fixed", fixedBreaks=den.break1, dataPrecision=5)
color.code.fitted.car = findColours(class.fitted.car, color.pallete)

class.fitted.sar = classIntervals(var=mn.poly$fitted.sar, n=6, style="fixed", fixedBreaks=den.break1, dataPrecision=5)
color.code.fitted.sar = findColours(class.fitted.sar, color.pallete)

leg.txt = c("0-5000", "5000-10000", "10000-20000","20000-40000","40000-80000",'80000-100000')

# plot
par(mfrow=c(2,1), oma = c(0,0,4,0) + 0.1, mar = c(0,0,1,0) + 0.1)

plot(mn.poly, col=color.code.fitted.car)
title("a)  CAR fitted rate")

plot(mn.poly, col=color.code.fitted.sar)
title("b)  SAR fitted rate")
legend("bottomright", legend=leg.txt, cex=0.8, bty="n", horiz = FALSE, fill = color.pallete)



