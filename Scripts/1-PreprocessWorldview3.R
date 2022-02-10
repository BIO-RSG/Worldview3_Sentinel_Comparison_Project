rm(list=ls())
library(raster)
library(XML)

#Import files
#this will only run on the multispectral files
#Raster data
ras.data = "C:\\Users\\wilsonkri\\Documents\\SatelliteData\\Worldview-Data-ESI\\012266432010_01\\012266432010_01_P001_MUL\\19AUG11151735-M2AS-012266432010_01_P001.tif"
#Metadata file
meta.file = "C:\\Users\\wilsonkri\\Documents\\SatelliteData\\Worldview-Data-ESI\\012266432010_01\\012266432010_01_P001_MUL\\19AUG11151735-M2AS-012266432010_01_P001.xml"
##Data to output
out.ras = "C:\\Users\\wilsonkri\\Documents\\SatelliteData\\Worldview-Data-ESI\\Processed\\20190811\\TOA-reflectance-20190811.tif"
png.out = "C:\\Users\\wilsonkri\\Documents\\SatelliteData\\Worldview-Data-ESI\\Processed\\20190811\\TOA-reflectance-20190811.png"
#

#Define sensor specific gain and offset (these may change with time)
#2018v0
WV3gain = c(0.938,#coastal blue
            0.946,#blue
            0.958,#green
            0.979,#yellow
            0.969,#red
            1.027,#red-edge
            0.977,#NIR1
            1.007)#NIR2
WV3offset = c(-13.099,#coastal blue
            -9.409,#blue
            -7.771,#green
            -5.489,#yellow
            -4.579,#red
            -5.552,#red-edge
            -6.508,#NIR1
            -3.699)#NIR2

##Define sensor specific Esun constant: shouldn't need to change these
#Thuiller 2003
Esun = c(1757.890,#coastal blue
         2004.610,#blue
         1830.180,#green
         1712.070,#yellow
         1535.330,#red
         1348.080,#red-edge
         1055.940,#NIR1
         858.770)#NR2
##

#Extract Data out of the metadata information
meta.file = xmlParse(meta.file)
meta.file = xmlToList(meta.file)
#Data for conversion to TOA Spectral Radiance
coastalblue = as.numeric(c(meta.file$IMD$BAND_C$ABSCALFACTOR, meta.file$IMD$BAND_C$EFFECTIVEBANDWIDTH))
blue = as.numeric(c(meta.file$IMD$BAND_B$ABSCALFACTOR, meta.file$IMD$BAND_B$EFFECTIVEBANDWIDTH))
green = as.numeric(c(meta.file$IMD$BAND_G$ABSCALFACTOR, meta.file$IMD$BAND_G$EFFECTIVEBANDWIDTH))
yellow = as.numeric(c(meta.file$IMD$BAND_Y$ABSCALFACTOR, meta.file$IMD$BAND_Y$EFFECTIVEBANDWIDTH))
red = as.numeric(c(meta.file$IMD$BAND_R$ABSCALFACTOR, meta.file$IMD$BAND_R$EFFECTIVEBANDWIDTH))
rededge = as.numeric(c(meta.file$IMD$BAND_RE$ABSCALFACTOR, meta.file$IMD$BAND_RE$EFFECTIVEBANDWIDTH))
nir1 = as.numeric(c(meta.file$IMD$BAND_N$ABSCALFACTOR, meta.file$IMD$BAND_N$EFFECTIVEBANDWIDTH))
nir2 = as.numeric(c(meta.file$IMD$BAND_N2$ABSCALFACTOR, meta.file$IMD$BAND_N2$EFFECTIVEBANDWIDTH))
radianc.conv = rbind(coastalblue,blue,green,yellow,red,rededge,nir1,nir2)
colnames(radianc.conv) = c("AbsScaleFactor","EffectiveBandwidth")
rm(coastalblue,blue,green,yellow,red,rededge,nir1,nir2)
#Data for conversion to TOA spectral reflectence
#Extract date info
julian.day = (meta.file$IMD$MAP_PROJECTED_PRODUCT$EARLIESTACQTIME)
out = strsplit(julian.day, "[T]")
out2 = strsplit(out[[1]][1], "[-]")
out.year = as.numeric(out2[[1]][1])
out.month = as.numeric(out2[[1]][2])
out.date = as.numeric(out2[[1]][3])
if (out.month<3){out.year=out.year-1
out.month = out.month+12}#If month is January or February it must be modified
out2 = strsplit(out[[1]][2], "[:]")
out.hour = as.numeric(out2[[1]][[1]])
out.min = as.numeric(out2[[1]][[2]])
out2 = strsplit(out2[[1]][3], "[Z]")
out.sec = as.numeric(out2[[1]][[1]])
#Calculate the Julian Day
UT = out.hour+(out.min/60.0)+(out.sec/3600.0)
A = trunc(out.year/100)
A2 = trunc(A/4)
B = 2-A+A2
JD = trunc(365.25*(out.year+4716))+trunc(30.6001*(out.month+1))+out.date+(UT/24)+B-1524.5
#Calculate Earth Sun Distance
D = JD-2451545.0
G = 357.529 + (0.98560028*D)
G = G*(pi/180)#convetr to radians
D.es = 1.00014 - (0.01671*cos(G)) - (0.00014*cos(2*G))#earth
if(D.es<0.983 || D.es>1.017){print("Failed Des calc: Value of of range")}
rm(julian.day,out,out2,out.year,out.month,out.date,out.hour,out.min,out.sec,UT,A,B,JD,D,G,A2)
#
cos.sun = 90-as.numeric(meta.file$IMD$IMAGE$MEANSUNEL)#solar zenith angle in degrees
cos.sun = cos.sun*(pi/180)#solar zenith angle in radians
cos.sun = cos(cos.sun)#final solar zenith angle to input
rm(meta.file)
#

#Read in Raster Data
ras.data = brick(ras.data)
names(ras.data) = c("cb","b","g","y", "r", "re","n1","n2")
ras.data = reclassify(ras.data, c(-Inf,0,NA))
ras.data = ras.data+min(minValue(ras.data))
for (i in 1:8){
  print(i)
  #convert to TOA radiance
  ras.data[[i]] = (WV3gain[i]*(ras.data[[i]]*(radianc.conv[i,"AbsScaleFactor"]/radianc.conv[i,"EffectiveBandwidth"]))) + WV3offset[i]
  #convert to TOA reflectance
  ras.data[[i]] = (ras.data[[i]]*(D.es^2)*pi)/(Esun[i]*cos.sun)
  }
ras.data

#Look at Data
png(png.out, width = 8, height=8, units = "in", res = 300)
list.col = c("blueviolet","blue3", "springgreen2","gold","orangered3", "red4", "grey", "black")
list.name = c("coastalblue", "blue","green","yellow","red","red edge","nir1","nir2")
par(mfrow = c(2,4),mar=c(4,4,3,1))
for (i in 1:8){density(ras.data[[i]],main=list.name[i],xlab="Reflectance",ylab="Density",col = list.col[i])}
dev.off()
#
writeRaster(ras.data,out.ras,format="GTiff",NAflag = NaN,overwrite=T)


