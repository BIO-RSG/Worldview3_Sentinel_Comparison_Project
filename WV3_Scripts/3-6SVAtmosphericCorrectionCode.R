##
rm(list=ls())
library(raster)
#
setwd("C:\\Users\\wilsonkri\\Documents\\Backup-Worldview-Data\\20190817\\")
##Import WV-2 data
sat.dat =brick("TOA-reflectance-20190817.tif")
names(sat.dat) = c("cb","b","g","y", "r", "re","n1","n2")
#
sixSV.dat = "./AC/acv2/ac097/output_"#link to 6SV output
#
out.ras = "BOA-Level2-reflectance-20190817.tif"
#
land.mask = raster("./ExtraFiles/LandMask.tif")
cloud.mask1 = shapefile("./ExtraFiles/cloud/CloudMask_Polygon.shp")
cloud.mask2 = shapefile("./ExtraFiles/cloud/cloud2_Polygon.shp")

##
sat.dat = mask(sat.dat,land.mask)
sat.dat = mask(sat.dat,cloud.mask1,inverse=T)
sat.dat = mask(sat.dat,cloud.mask2,inverse=T)
rm(land.mask,cloud.mask1,cloud.mask2)
#min.toa = minValue(sat.dat)
#names(min.toa) = c("cb","b","g","y", "r", "re","n1","n2")
#
##Loop through
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
par(mfrow = c(2,4),mar=c(3,3,1,1),mgp=c(2,1,0))
for(i in 1:8){hist(sat.dat[[i]],breaks=c(-0.1,seq(0,0.58,0.01)),ylim=c(0,2))}

plotRGB(sat.dat,r=5,g=3,b=2,ext=c(530000,534000,4958000,4960000),stretch="lin")
plot(mask.bright,add=T,alpha=0.5)

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
                                               0.1,1,0),byrow=T,nrow=2),right=F)
mask.n2 = reclassify(sat.dat$n2, matrix(c(0,0.2,1,
                                          0.1,1,0),byrow=T,nrow=2),right=F)
mask.bright = mask.blue*mask.red*mask.yellow*mask.green*mask.cb*mask.rededge*mask.n2#*mask.n1
rm(mask.blue,mask.red,mask.yellow,mask.green,mask.cb,mask.n2,mask.rededge)#,mask.n1)
#
mask.bright = reclassify(mask.bright,matrix(c(-0.5,0.5,NA,0.5,1.5,1),byrow=T,nrow=2))
sat.dat = mask(sat.dat, mask.bright)
sat.dat
#
writeRaster(sat.dat,out.ras,format="GTiff",NAflag = NaN,overwrite=T)
#

#getneg = function(x) {ifelse(any(x<0)==T,NA,1)}
#rasterOptions(maxmemory=7e+09,chunksize = 4e+08)
#beginCluster()  
#neg.mask = clusterR(sat.dat, calc, args=list(fun=getneg))#this takes awhile to run  
#endCluster() 
#plot(neg.mask,col=c("red"))
#sat.dat= mask(x=sat.dat, mask = neg.mask)
##

# Data Distribution -------------------------------------------------------

toa = sat.dat
boa = sat.dat
#

list.col = c("blueviolet","blue3", "springgreen2","gold","orangered3", "red4", "grey", "black")
list.name = c("coastalblue", "blue","green","yellow","red","red edge","nir1","nir2")
list.max = c(0.15,0.15,0.15,0.1,0.1,0.1,0.05,0.05)

#Look at Data
png("./AC/acv2/Compare-TOA-BOA-017.png", width = 8, height=6, units = "in", res = 300)
par(mfrow = c(2,4),mar=c(3,3,1,1),mgp=c(2,1,0))
for (i in 1:8){
  print(i)
  toa.den = density(toa[[i]],plot=F)
  boa.den = density(boa[[i]],plot=F)
  plot(-Inf, xlim=c(-0.01,list.max[i]),ylim=c(0,max(toa.den$y,boa.den$y)+5), xaxs = "i",yaxs="i", xlab="Reflectance",ylab="Density")
  abline(v=0)
  lines(toa.den,col = list.col[i],lty=2)
  lines(boa.den,col = list.col[i])
  legend("topright", legend=c("TOA","BOA"), col = list.col[i], lty=c(2,1),bty="n",title=list.name[i],seg.len =1,cex=0.8)
  }
dev.off()
#

#min.boa = minValue(sat.dat)
#names(min.boa) = c("cb","b","g","y", "r", "re","n1","n2")
png("./AC/acv2/minPixel.png",width=6,height=4,res=300,units="in")
par(mfrow=c(2,3),mar=c(3,3,1,1),mgp=c(2,1,0))
plot(min.toa, ylim=c(-0.01,0.17),pch=20, xlab="Band",ylab="Min Reflectance",main="AOT=0.017 (N1)")
points(min.boa.017,pch=20,col=rgb(1,0,0,0.5))
abline(h=0,lty=2)
plot(min.toa, ylim=c(-0.01,0.17),pch=20, xlab="Band",ylab="Min Reflectance",main="AOT=0.091 (N2)")
points(min.boa.091,pch=20,col=rgb(0,1,0,0.5))
abline(h=0,lty=2)
plot(min.toa, ylim=c(-0.01,0.17),pch=20, xlab="Band",ylab="Min Reflectance",main="AOT=0.092 (N2)")
points(min.boa.092,pch=20,col=rgb(0,0,1,0.5))
abline(h=0,lty=2)
plot(min.toa, ylim=c(-0.01,0.17),pch=20, xlab="Band",ylab="Min Reflectance",main="AOT=0.097 (RE)")
points(min.boa.097,pch=20,col=rgb(1,0,0,0.5))
abline(h=0,lty=2)
plot(min.toa, ylim=c(-0.01,0.17),pch=20, xlab="Band",ylab="Min Reflectance",main="AOT=0.124 (R)")
points(min.boa.124,pch=20,col=rgb(0,1,0,0.5))
abline(h=0,lty=2)
plot(min.toa, ylim=c(-0.01,0.17),pch=20, xlab="Band",ylab="Min Reflectance",main="AOT=0.128 (R)")
points(min.boa.128,pch=20,col=rgb(0,0,1,0.5))
abline(h=0,lty=2)
dev.off()
#
