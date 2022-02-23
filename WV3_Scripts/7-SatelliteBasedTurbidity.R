library(raster)
library(ncdf4)
library(remotes)
#install_github("BIO-RSG/oceancolouR")
library(oceancolouR)
library(pals)
library(palr)
options(scipen = 999)

##Read in your data
setwd("C:\\Where\\is\\my\\data\\")
wv.dat = brick("./StripeCorrection/PD3-use/BOA-destripe.tif")
names(wv.dat) = c("cb","b","g","y","r","re","n1","n2")
depth = raster("./ExtraFiles/Depth.tif")
depth[depth<10]=NA
out.folder = c("./WaterQualityDS//")
##

###chl oc
in.dat = stack(wv.dat$g,wv.dat$cb)
names(in.dat) = c("Rrs_561","Rrs_443")
chl.oc2 = ocx(in.dat,
              c("Rrs_443"),
              c("Rrs_561"),
              c(0.1977, -1.8117,1.9743, -2.5635, -0.7218),#Landsat, https://oceancolor.gsfc.nasa.gov/atbd/chlor_a/
              use_443nm=T)
chl.oc2[chl.oc2<=0]=NA
chl.oc2[chl.oc2>1000]=NA
chl.oc2.m = mask(chl.oc2,depth)
writeRaster(chl.oc2, paste0(out.folder,"chl-oc2.tif"),format="GTiff",NAflag = NaN,overwrite=T)
png(paste0(out.folder,"chl-oc2.png"),res=300,width=8,height=3,units = "in")
par(mfcol=c(1,3),mar=c(3,3,0,0),oma=c(0,0,0.5,2),mgp=c(2,1,0))
density(chl.oc2.m,xlab="chla (ug/l)",log="x",col=rgb(0,0,1,0.5),lwd=2,xlim=c(minValue(chl.oc2),maxValue(chl.oc2)))
a=density(chl.oc2,plot=F)
lines(a,col=rgb(1,0,0,0.5),lwd=2)
legend("topright",legend=c("All Points", ">10m Points"),col=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),lty=1,lwd=2,bty="n")
boxplot(na.omit(getValues(chl.oc2)),ylab="chla (ug/l)",xlab="oc2",log="y",
        ylim=c(minValue(chl.oc2),maxValue(chl.oc2)),xlim=c(0.75,1.75),at=1,col=rgb(1,0,0,0.5))
boxplot(na.omit(getValues(chl.oc2.m)),log="y",add=T,at=1.5,col=rgb(0,0,1,0.5))
plot(chl.oc2, breaks = c(0,seq(0.05,0.4,0.01),seq(0.5,1,0.1),ceiling(maxValue(chl.oc2))),
     col= chl_pal (44), legend=F,maxpixels = ncell(chl.oc2))
plot(chl.oc2, breaks = c(0,seq(0.05,0.4,0.01),seq(0.5,1,0.1),ceiling(maxValue(chl.oc2))),
     col= chl_pal (44),legend.only=T, smallplot=c(0.9,0.94, 0.05,0.95),
     horizontal=F,legend.args=list(text="chla",line=0.5,cex=0.7))
dev.off()
rm(in.dat,a,chl.oc2.m,chl.oc2)
##

###chl_oc3
in.dat = stack(wv.dat$g,wv.dat$b,wv.dat$cb)
names(in.dat) = c("Rrs_561","Rrs_482","Rrs_443")
chl.oc3 = ocx(in.dat,
              c("Rrs_443","Rrs_482"),
              c("Rrs_561"),
              c( 0.2412, -2.0546, 1.1776, -0.5538, -0.4570),
              use_443nm=T)
chl.oc3[chl.oc3<=0]=NA
chl.oc3[chl.oc3>1000]=NA
chl.oc3.m = mask(chl.oc3,depth)
writeRaster(chl.oc3, paste0(out.folder,"chl-oc3.tif"),format="GTiff",NAflag = NaN,overwrite=T)
png(paste0(out.folder,"chl-oc3.png"),res=300,width=8,height=3,units = "in")
par(mfcol=c(1,3),mar=c(3,3,0,0),oma=c(0,0,0.5,2),mgp=c(2,1,0))
density(chl.oc3.m,xlab="chla (ug/l)",log="x",col=rgb(0,0,1,0.5),lwd=2,xlim=c(minValue(chl.oc3),maxValue(chl.oc3)))
a=density(chl.oc3,plot=F)
lines(a,col=rgb(1,0,0,0.5),lwd=2)
legend("topright",legend=c("All Points", ">10m Points"),col=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),lty=1,lwd=2,bty="n")
boxplot(na.omit(getValues(chl.oc3)),ylab="chla (ug/l)",xlab="oc3",log="y",
        ylim=c(minValue(chl.oc3),maxValue(chl.oc3)),xlim=c(0.75,1.75),at=1,col=rgb(1,0,0,0.5))
boxplot(na.omit(getValues(chl.oc3.m)),log="y",add=T,at=1.5,col=rgb(0,0,1,0.5))
plot(chl.oc3, breaks = c(0,seq(0.05,0.4,0.01),seq(0.5,1,0.1),ceiling(maxValue(chl.oc3))),
     col= chl_pal (44), legend=F,maxpixels = ncell(chl.oc3))
plot(chl.oc3, breaks = c(0,seq(0.05,0.4,0.01),seq(0.5,1,0.1),ceiling(maxValue(chl.oc3))),
     col= chl_pal (44),legend.only=T, smallplot=c(0.9,0.94, 0.05,0.95),
     horizontal=F,legend.args=list(text="chla",line=0.5,cex=0.7))
dev.off()
rm(in.dat,a,chl.oc3.m,chl.oc3)
##

##chl re gons
#define coefficients
#astar_chl=0.015	
#a0=1.61
#a1=0.082
#a2=0.6
#a3=0.7
#a4 = 0.40
#a5= 1.05	
#val1 = 0.005
#val2 = 0.63	
#red.664 = wv.dat$r #centre wavelength 660
#rededge.704 = wv.dat$re #centre wavelength 722
# this band is bad in general
#nir.782 = wv.dat$n1 #centre wavelength 824
#
## Calculate
#bb = (a0*nir.782)/(a1-a2*nir.782)
#rm = rededge.704/red.664
#chl.gons = (rm*(a3+bb)) - a4 - bb^a5
#chl.gons = chl.gons/astar_chl
#chl.gons[chl.gons<0]=NA #all negative to NA
#val.mask = reclassify(red.664, matrix(c(-Inf,val1,NA, val1,Inf,1),byrow=T,nrow=2),right=F)#set validity range on input bands
#val.mask2 = reclassify(nir.782/rededge.704, matrix(c(-Inf,val2,NA, val2,Inf,1),byrow=T,nrow=2),right=F)#set validity range on input bands
#val.mask = val.mask*val.mask2
#chl.gons = mask(x=chl.gons, mask = val.mask)
#writeRaster(chl.gons, paste0(out.folder,"chl-gons.tif"),format="GTiff",NAflag = NaN,overwrite=T)
#rm(red.664,rededge.704,nir.782,bb,rm,a0,a1,a2,a3,a4,a5,val1,val2,astar_chl,val.mask,val.mask2)

##chl re moses
#define coefficients
#a1 = 232.29
#a2 = 23.173
#red.664 = wv.dat$r #centre wavelength 660
#rededge.704 = wv.dat$re #centre wavelength 722
# this band is bad in general
#nir.780 = wv.dat$n1 #centre wavelength 824
#calculate
#chl.moses = a1*( (red.664^-1 - rededge.704^-1) * nir.780 ) + a2
#chl.moses[chl.moses<0]=NA
#writeRaster(chl.moses, paste0(out.folder,"chl-moses.tif"),format="GTiff",NAflag = NaN,overwrite=T)
#rm(a1,a2,red.664,rededge.704,nir.780)
##

#SPM nechad
#WV-3 values values
#Nechad et al 2010 "Calibration and validation of a generic multisensor algorithm for mapping of totalsuspended matter in turbid waters"
A = 327.84 #centre 661, use 660, Table 4
C = 0.1708 #centre 661, use 660, Table 1
red.664 = wv.dat$r #centre wavelength 661
#calculate
spm.nechad = (A * red.664)/(1.-red.664/C)
##
spm.nechad[spm.nechad<=0]=NA
spm.nechad[spm.nechad>1000]=NA
spm.nechad.m = mask(spm.nechad,depth)
writeRaster(spm.nechad, paste0(out.folder,"spm-nechad.tif"),format="GTiff",NAflag = NaN,overwrite=T)
png(paste0(out.folder,"spm-nechad.png"),res=300,width=8,height=3,units = "in")
par(mfcol=c(1,3),mar=c(3,3,0,0),oma=c(0,0,0.5,2),mgp=c(2,1,0))
density(spm.nechad.m,xlab="spm (gm-3)",log="x",col=rgb(0,0,1,0.5),lwd=2,xlim=c(minValue(spm.nechad),maxValue(spm.nechad)))
a=density(spm.nechad,plot=F)
lines(a,col=rgb(1,0,0,0.5),lwd=2)
legend("topright",legend=c("All Points", ">10m Points"),col=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),lty=1,lwd=2,bty="n")
boxplot(na.omit(getValues(spm.nechad)),ylab="spm (gm-3)",xlab="Nechad",log="y",
        ylim=c(minValue(spm.nechad),maxValue(spm.nechad)),xlim=c(0.75,1.75),at=1,col=rgb(1,0,0,0.5))
boxplot(na.omit(getValues(spm.nechad.m)),log="y",add=T,at=1.5,col=rgb(0,0,1,0.5))
plot(spm.nechad, breaks = c(0,seq(1,5,0.1),seq(6,10,1),ceiling(maxValue(spm.nechad))),
     col=ocean.deep (49), legend=F,maxpixels = ncell(spm.nechad))
plot(spm.nechad, breaks = c(0,seq(1,5,0.1),seq(6,10,1),ceiling(maxValue(spm.nechad))),
     col= ocean.deep (49),legend.only=T, smallplot=c(0.9,0.94, 0.05,0.95),
     horizontal=F,legend.args=list(text="spm",line=0.5,cex=0.7))
dev.off()
rm(A,C,red.664,a,spm.nechad,spm.nechad.m)
##

#turbidity dogliotti
##values for wv-3 wavelengths based on the centre wavelengths
#WV-3 centre wavelength in order red, red-edge, nir1: 661,724,832
#Nechad et al 2009 Calibration and validation of a generic multisensor algorithm for mapping of turbidity in coastal waters
a.t.red = 261.11 #centre 661, use 660
#If use red edge band
a.t.nir =  679.38 #centre 724 use 725
#If use nir1 band
#a.t.nir =  1613.09 #centre 832 use 832.5
#Nechad et al 2010 "Calibration and validation of a generic multisensor algorithm for mapping of totalsuspended matter in turbid waters"
c.t.red = 0.1708 #centre 661, use 660
#If use red edge band
c.t.nir = 0.1937 #centre 724 use 725
#If use nir1 band
#c.t.nir = 0.2101 #centre 832 use 832.5
##
low.lim = 0.05
up.lim =  0.07
red.664 = wv.dat$r #centre wavelength 660
# this band is bad in general
nir.833 = wv.dat$n1 #centre wavelength 824
#find index of most turbid with nir band
red.index = reclassify(red.664,right=F, matrix(c(-Inf,up.lim,0,up.lim,Inf,1),byrow=T,nrow=2))
#find index of regions to blend in between
blend.index = reclassify(red.664,right=F, matrix(c(-Inf,low.lim,0,
                                  low.lim,up.lim,1,up.lim,Inf,0),byrow=T, nrow=3))
#calculate turbidity
turb.red = (a.t.red * red.664)/(1.0 - red.664/c.t.red)
turb.nir = (a.t.nir*nir.833)/(1.0 - nir.833/c.t.nir)
w = (red.664 - low.lim)/(up.lim - low.lim)
turb.blend = ((1.0-w) *turb.red) + (w * turb.nir)
#Fill final layer
turb.d = turb.red#define a raster of only 0
turb.d[red.index==1]=turb.nir[red.index==1]#replace high values
turb.d[blend.index==1]=turb.blend[blend.index==1]#replace high values
turb.d[turb.d<0]=NA
turb.d[turb.d>1000]=NA
turb.d.m = mask(turb.d,depth)
writeRaster(turb.d , paste0(out.folder,"turb-dogliotti.tif"),format="GTiff",NAflag = NaN,overwrite=T)
#
png(paste0(out.folder,"turb-dogliotti.png"),res=300,width=8,height=3,units = "in")
par(mfcol=c(1,3),mar=c(3,3,0,0),oma=c(0,0,0.5,2),mgp=c(2,1,0))
density(turb.d.m,xlab="FNU",log="x",col=rgb(0,0,1,0.5),lwd=2,xlim=c(minValue(turb.d),maxValue(turb.d)),ylim=c(0,8))
a=density(turb.d,plot=F)
lines(a,col=rgb(1,0,0,0.5),lwd=2)
legend("topright",legend=c("All Points", ">10m Points"),col=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),lty=1,lwd=2,bty="n")
boxplot(na.omit(getValues(turb.d)),ylab="FNU",xlab="Turbidity",log="y",
        ylim=c(minValue(turb.d),maxValue(turb.d)),xlim=c(0.75,1.75),at=1,col=rgb(1,0,0,0.5))
boxplot(na.omit(getValues(turb.d.m)),log="y",add=T,at=1.5,col=rgb(0,0,1,0.5))
plot(turb.d, breaks = c(0,seq(1,5,0.1),seq(6,10,1),ceiling(maxValue(turb.d))),
     col=ocean.deep (49), legend=F,maxpixels = ncell(turb.d))
plot(turb.d, breaks = c(0,seq(1,5,0.1),seq(6,10,1),ceiling(maxValue(turb.d))),
     col= ocean.deep (49),legend.only=T, smallplot=c(0.9,0.94, 0.05,0.95),
     horizontal=F,legend.args=list(text="FNU",line=0.5,cex=0.7))
dev.off()
###
rm(a.t.red,c.t.red,a.t.nir,c.t.nir,low.lim,up.lim,red.664,nir.833)
rm(red.index,blend.index,turb.red,turb.nir,w,turb.blend,a)
rm(turb.d,turb.d.m)
##




