##
#This script will atmospherically correct the satellite data
#And mask your data
library(raster)
#
setwd("C:\\Where\\is\\my\\data\\")
##Import WV-2 data
sat.dat =brick("TOA-reflectance-20190817.tif")
names(sat.dat) = c("cb","b","g","y", "r", "re","n1","n2")
#
sixSV.dat = "./OutData/output_"#link to 6SV output created in 2b
#
out.ras = "BOA-Level2-reflectance-20190817.tif"
#
land.mask = raster("./ExtraFiles/LandMask.tif")#Predefined
cloud.mask1 = shapefile("./ExtraFiles/cloud/CloudMask_Polygon.shp")#Predefined
cloud.mask2 = shapefile("./ExtraFiles/cloud/cloud2_Polygon.shp")#Predefined 

##Land and cloud masking
sat.dat = mask(sat.dat,land.mask)
sat.dat = mask(sat.dat,cloud.mask1,inverse=T)
sat.dat = mask(sat.dat,cloud.mask2,inverse=T)
rm(land.mask,cloud.mask1,cloud.mask2)

#
#Atmospheric correction
for (i in 1:dim(sat.dat)[3]){
  print(i)
  sixSV.dat.in = paste(sixSV.dat,i,".txt",sep="")
  sixSV.dat.in = read.delim(sixSV.dat.in,header=F,colClasses = "character",fill=T,sep=c(""),skipNul = T,na.strings = "")
  a.ref = which(sixSV.dat.in == "xap", arr.ind=TRUE)#find position of ref coefficients 
  a.ref=sixSV.dat.in[(a.ref[1]+1),]#extract row of data
  a.ref = a.ref[, !apply(is.na(a.ref), 2, all)]
 #Generate coeffificients 
  if (a.ref[1]!=":"){
    xap = as.numeric(a.ref[1])
    xb = as.numeric(a.ref[2])
    xc = as.numeric(a.ref[3])}
  if (a.ref[1]==":"){
    xap = as.numeric(a.ref[2])
    xb = as.numeric(a.ref[3])
    xc = as.numeric(a.ref[4])}
 #Calculate BOA reflectance
  y.ref=(xap*sat.dat[[i]])-xb #6SV formula
  sat.dat[[i]]=y.ref/(1.+(xc*y.ref))#6SV formula
}
rm(y.ref,a.ref,sixSV.dat,sixSV.dat.in,i,xap,xb,xc)
#
##
#Generate bright pixel mask
mask.cb= reclassify(sat.dat$cb, matrix(c(0,0.1,1,
                                         0.1,1,0),byrow=T,nrow=2),right=F)
mask.blue = reclassify(sat.dat$b, matrix(c(0,0.1,1,
                                           0.1,1,0),byrow=T,nrow=2),right=F)
mask.green = reclassify(sat.dat$g, matrix(c(0,0.1,1,
                                            0.1,1,0),byrow=T,nrow=2),right=F)
mask.yellow = reclassify(sat.dat$y, matrix(c(0,0.1,1,
                                             0.1,1,0),byrow=T,nrow=2),right=F)
mask.red = reclassify(sat.dat$r, matrix(c(0,0.1,1
                                          ,0.1,1,0),byrow=T,nrow=2),right=F)
mask.rededge = reclassify(sat.dat$re, matrix(c(0,0.2,1,
                                               0.2,1,0),byrow=T,nrow=2),right=F)
mask.n2 = reclassify(sat.dat$n2, matrix(c(0,0.2,1,
                                          0.2,1,0),byrow=T,nrow=2),right=F)
mask.bright = mask.blue*mask.red*mask.yellow*mask.green*mask.cb*mask.rededge*mask.n2#*mask.n1
rm(mask.blue,mask.red,mask.yellow,mask.green,mask.cb,mask.n2,mask.rededge)#,mask.n1)
#
mask.bright = reclassify(mask.bright,matrix(c(-0.5,0.5,NA,0.5,1.5,1),byrow=T,nrow=2))
sat.dat = mask(sat.dat, mask.bright)
sat.dat
#
#Final atmospherically corrected and masked data
writeRaster(sat.dat,out.ras,format="GTiff",NAflag = NaN,overwrite=T)
#




