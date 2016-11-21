# GSM spatial
# Data Merge and Clean

setwd('C:/Users/Dannie/Dropbox/Dissertation/GMSData/R/SpatialAnalysis')

library(geoR)
library(akima)
library(gstat)

#################################
# merge five hour data into one file 
{
  D.raw<-read.table('fishnetpoint05.txt',sep=',',header = T)
  D05<-D.raw[c('long','lat','freq')]# other hour
  D05$hour<-rep(5,nrow(D05))
  
  D.raw<-read.table('fishnetpoint06.txt',sep=',',header = T)
  D06<-D.raw[c('long','lat','freq')]# other hour
  D06$hour<-rep(6,nrow(D06))
  
  D.raw<-read.table('fishnetpoint07.txt',sep=',',header = T)
  D07<-D.raw[c('long','lat','freq')]# other hour
  D07$hour<-rep(7,nrow(D07))
  
  D.raw<-read.table('fishnetpoint08.txt',sep=',',header = T)
  D08<-D.raw[c('long','lat','freq')]# other hour
  D08$hour<-rep(8,nrow(D08))
  
  D.raw<-read.table('fishnetpoint09.txt',sep=',',header = T)
  D09<-D.raw[c('long','lat','freq')]# other hour
  D09$hour<-rep(9,nrow(D09))
  
  D<-rbind(D05,D06,D07,D08,D09)
  
  D05<-NULL
  D06<-NULL
  D07<-NULL
  D08<-NULL
  D09<-NULL
  D.raw<-NULL
  
  write.csv(D,'GSM_all.csv')
}
# end of merge data sets into one file


# 

filelist<-list.files(pattern = '*.txt')
filelist<-filelist[seq(2,15,3)]

D<-data.frame()
for (i in 1:5)
{
  D.raw<-read.table(filelist[i],sep=',',header = T)
  D05<-D.raw[c('FID','long','lat','freq','D2Metro','D2Road')]# other hour
  D05$hour<-rep(i+4,nrow(D05))
  D<-rbind(D,D05)
  rm(D05)
  rm(D.raw)
}

write.csv(D,'GSM_all_withDist.csv')
