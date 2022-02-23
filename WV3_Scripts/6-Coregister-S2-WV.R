library(raster)
library(RStoolbox)


wv.dat = brick("C:\\Where\\is\\my\\data//BOA-reflectance-20190811.tif")
names(wv.dat) = c("cb","b","g","y","r","re","n1","n2")
s2.dat = stack("C:\\Where\\is\\my\\data\\20160913Level2A.tif")
names(s2.dat) = c("b","g","r","n","depth","ndvi","gndvi","class")

#Crop the S-2 image to the WV image for speed
s2.dat=crop(s2.dat,extent(wv.dat))
wv.dat=crop(wv.dat,extent(s2.dat))

#Change the S-2 data to have 2-m pixel resolution
s2.dat = disaggregate(s2.dat, fact=5)
s2.dat=crop(s2.dat,extent(wv.dat))
wv.dat=crop(wv.dat,extent(s2.dat))

#Reproject the S-2 data to have the same projection, don't use project raster messes with #
s2.blue = wv.dat[[1]]
s2.blue[] =-1 
s2.blue[] = as.matrix(s2.dat$b)
s2.blue1 = wv.dat[[1]]
s2.blue1[] =-1 
s2.blue1[] = as.matrix(s2.dat$g)
s2.blue2 = wv.dat[[1]]
s2.blue2[] =-1 
s2.blue2[] = as.matrix(s2.dat$r)
s2.blue3 = wv.dat[[1]]
s2.blue3[] =-1 
s2.blue3[] = as.matrix(s2.dat$depth)
s2.blue5 = wv.dat[[1]]
s2.blue5[] =-1 
s2.blue5[] = as.matrix(s2.dat$ndvi)
s2.blue4 = wv.dat[[1]]
s2.blue4[] =-1 
s2.blue4[] = as.matrix(s2.dat$class)
s2.dat = stack(s2.blue,s2.blue1,s2.blue2,s2.blue3,s2.blue4,s2.blue5)
names(s2.dat) = c("b","g","r","depth","class","ndvi")
rm(s2.blue,s2.blue1,s2.blue2,s2.blue3,s2.blue4)

#write out in batches and delete later
writeRaster(s2.dat[[1:2]],
            "C:\\Where\\is\\my\\data//S-2Comparison//S2-v2//s2-2m-reprj-tmp1.tif",
            format="GTiff",NAflag = NaN,overwrite=T)
writeRaster(s2.dat[[3:5]],
            "C:\\Where\\is\\my\\data\\S-2Comparison//S2-v2//s2-2m-reprj-tmp2.tif",
            format="GTiff",NAflag = NaN,overwrite=T)
s2.dat = stack("C:\\Where\\is\\my\\data\\S-2Comparison//S2-v2//s2-2m-reprj-tmp1.tif",
               "C:\\Where\\is\\my\\data\\S-2Comparison//S2-v2//s2-2m-reprj-tmp2.tif")

#Write out only downscaled and repojected data
writeRaster(s2.dat,
            "C:\\Where\\is\\my\\data//S-2Comparison//S2-v2//s2-2m-reprj.tif",
            format="GTiff",NAflag = NaN,overwrite=T)
s2.dat =brick( "C:\\Where\\is\\my\\data\\S-2Comparison//S2-v2//s2-2m-reprj.tif")
names(s2.dat) = c("b","g","r","depth","class")

plot(s2.dat$b,xlim=c(527400,527420),ylim=c(4956080,4956100))
plot(wv.dat$b,xlim=c(527400,527420),ylim=c(4956080,4956100),add=T,alpha=0.5)

#Shift the S-2 image to the worldview image
s2.shift = coregisterImages(slave=s2.dat[[c("b","g","r")]],
                            master=wv.dat[[c("b","g","r")]],
                            shift=5,
                            shiftInc = 0.5, #Times by resolution of raster so to go up by 1 need to set to 0.5
                            nSamples = 1e+05, reportStats = T)
write.csv(s2.shift$MI, row.names = F,
          "C:\\Where\\is\\my\\data\\S-2Comparison\\BestS2Shift-Destripe-BGR-fixprj.txt")

#Does the destripe versus uncorrected suggest different shifts
library(pals)
ds.dat = read.csv("C:\\Where\\is\\my\\data\\S-2Comparison\\BestS2Shift-Destripe-BGR-fixprj.txt")
new.index = sort(ds.dat$mi, index=T,decreasing = T)
ds.dat=ds.dat[new.index[[2]],]
#boa.dat = read.csv("C:\\Users\\wilsonkri\\Documents\\Backup-Worldview-Data\\20190811\\S-2Comparison\\BestS2Shift-boa-BGR.txt")
#new.index = sort(boa.dat$mi, index=T,decreasing = T)
#boa.dat=boa.dat[new.index[[2]],]
par(mfcol=c(1,1),mar=c(3,3,1,1),mgp=c(2,1,0))
in.col =cubicl(n=length(unique(ds.dat$y)))
#plot(y=boa.dat$mi,x=boa.dat$x, pch=ifelse(boa.dat$y==0&boa.dat$x==0,5,20),
 #    ylim=c(0.74,0.82),main="BOA Uncorrected",col=in.col[as.numeric(as.factor(boa.dat$y))],
  #   xlab="X Shift", ylab="Mutual Information")
#legend("bottom", legend = unique(boa.dat$y),col=in.col,bty="n",pch=19,ncol=7,title="Y Shift",
#       x.intersp=0.75, text.width=1.5)
#legend("topleft", legend="x = 0 y = 0",pch=5,bty="n")
plot(y=ds.dat$mi,x=ds.dat$x,  pch=ifelse(ds.dat$y==0&ds.dat$x==0,5,20),
     ylim=c(0.76,0.84),main="BOA Corrected",col=in.col[as.numeric(as.factor(ds.dat$y))],
     xlab="X Shift", ylab="Mutual Information")
legend("top", legend = unique(ds.dat$y),col=in.col,bty="n",pch=19,ncol=7,title="Y Shift",
       x.intersp=0.75, text.width=1.5)
legend("bottom", legend="x = 0 y = 0",pch=5,bty="n")

#Zoom in on highest points
library(plotrix)
par(mfcol=c(1,1),mar=c(3,3,1,1),mgp=c(2,1,0))
#plot(y=boa.dat$mi,x=boa.dat$x, pch=ifelse(boa.dat$y==0&boa.dat$x==0,5,20),
 #    xlim=c(1,4),
  #   ylim=c(0.812,0.815),main="BOA Uncorrected",col=in.col[as.numeric(as.factor(boa.dat$y))],
   #  xlab="X Shift", ylab="Mutual Information")
#boa.dat$mi =round(boa.dat$mi,6)
#addtable2plot(x=1,y=0.814, table =boa.dat[1:5,],title="Best Shift")
plot(y=ds.dat$mi,x=ds.dat$x,  pch=ifelse(ds.dat$y==0&ds.dat$x==0,5,20),
     xlim=c(1,4),
     ylim=c(0.805,0.817),main="BOA Corrected",col=in.col[as.numeric(as.factor(ds.dat$y))],
     xlab="X Shift", ylab="Mutual Information")
ds.dat$mi =round(ds.dat$mi,6)
addtable2plot(x=1,y=0.814, table =ds.dat[1:5,],title="Best Shift")


#Coregister based on minimum shift with maximum MI for boa uncorrected wv data 
head(boa.dat)
#x=2 y=-4 is the smallest shift with maximum MI uncorrected
shift.x = 2
shift.y = -4
shift.matrix = matrix(c(shift.x,shift.y),ncol=2)/2
colnames(shift.matrix)= c("x","y")
#generate the raster dims
s2.shift = coregisterImages(slave=s2.dat[[c("b","g","r")]],
                            master=wv.dat[[c("b","g","r")]],
                            shift=shift.matrix,
                            #shiftInc = 0.5, #argument ignored
                            nSamples = 1e+05, reportStats = T)
#
s2.class.out = s2.shift$coregImg[[1]]
s2.class.out[] = as.matrix(s2.class)
writeRaster(stack(s2.shift$coregImg,s2.class.out),
            "C:\\Where\\is\\my\\data\\S-2Comparison//NewRF-05052020-2ndvegadd//s2-2m-reprj-shiftboa.tif",
            format="GTiff",NAflag = NaN,overwrite=T)

#Coregister based on minimum shift with maximum MI for boa destriped wv data
head(ds.dat)
#x=2 y=-4 is the smallest shift with maximum MI
shift.x = 2
shift.y = -4
shift.matrix = matrix(c(shift.x,shift.y),ncol=2)/2
colnames(shift.matrix)= c("x","y")
#generate the raster dims
s2.shift = coregisterImages(slave=s2.dat[[c("b","g","r")]],
                            master=wv.dat[[c("b","g","r")]],
                            shift=shift.matrix,
                            #shiftInc = 0.5, #argument ignored
                            nSamples = 1e+05, reportStats = T)
#
s2.class.out = s2.shift$coregImg[[1]]
s2.class.out[] = as.matrix(s2.dat$class)
s2.depth.out = s2.shift$coregImg[[1]]
s2.depth.out[] = as.matrix(s2.dat$depth)
s2.class.out[] = as.matrix(s2.dat$class)
s2.ndvi.out = s2.shift$coregImg[[1]]
s2.ndvi.out[] = as.matrix(s2.dat$ndvi)
shift.out = stack(s2.shift$coregImg,s2.depth.out,s2.class.out,s2.ndvi.out)
#Match extent
shift.out = crop(y=wv.dat,x=shift.out)
#Mask where WV data is out of bounds
wv.mask.ext = raster("C:\\Where\\is\\my\\data//Rasters//TOA-reflectance-20190811.tif")
wv.mask.ext[wv.mask.ext>0] = 1
wv.mask.ext = crop(y=shift.out,x=wv.mask.ext)
shift.out = mask(x=shift.out, mask=wv.mask.ext)

writeRaster(shift.out,
            "C:\\Where\\is\\my\\data\\S-2Comparison//S2-v2//s2-2m-reprj-shiftds.tif",
            format="GTiff",NAflag = NaN,overwrite=T)


