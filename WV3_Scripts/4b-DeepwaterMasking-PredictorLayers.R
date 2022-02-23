library(raster)
library(RStoolbox)

setwd("C:\\Where\\is\\my\\data\\")

wv.dat = brick("./Rasters/BOA-destripe.tif")
names(wv.dat) = c("cb","b","g","y","r","re","n1","n2")

#Generate 10m mask
depth.dat = raster("./Extrafiles/Depth.tif")
depth.dat[depth.dat>=10] = NA
writeRaster(depth.dat, "./Extrafiles/Depth10m-mask.tif", format="GTiff", NAflag = NaN)
depth.dat = raster("./Extrafiles/Depth10m-mask.tif")
###

##Mask data
wv.dat = mask(wv.dat,depth.dat)
writeRaster(wv.dat, "./Rasters/model-in/BOA-destripe-10mmask.tif", format="GTiff", NAflag = NaN)
wv.dat =brick("./Rasters/model-in/BOA-destripe-10mmask.tif")
names(wv.dat) = c("cb","b","g","y","r","re","n1","n2")

#Generate raster pca on 10m depth for the six bands
pca.dat = rasterPCA(wv.dat[[1:6]], spca=F)
writeRaster(pca.dat$map, "./Rasters/model-in/PCA-nostand-BOA-ds-6b-10mm.tif", format="GTiff", NAflag = NaN)
pca.dat =brick("./Rasters/model-in/PCA-nostand-BOA-ds-6b-10mm.tif")
names(pca.dat)=c("PC1","PC2","PC3","PC4","PC5","PC6")

##Import wcc
dii = brick("./NewWCC/Lyzenga1985-ds.tif")
names(dii) = c("diicbb","diicbg","diicby","diicbr", "diibg", "diiby","diibr", "diigy","diigr","diiyr")
dii = mask(dii,depth.dat)
writeRaster(dii, "./Rasters/model-in/Lyzenga1985-ds-10mmask.tif", format="GTiff", NAflag = NaN)
dii = brick("./Rasters/model-in/Lyzenga1985-ds-10mmask.tif")
names(dii) = c("diicbb","diicbg","diicby","diicbr", "diibg", "diiby","diibr", "diigy","diigr","diiyr")

##Import SDB
sdb = raster("./ExtraFiles/Depth/subset_BOA-destripe-bg_empBathymetry-mf.tif")
sdb  = mask(sdb,depth.dat)
writeRaster(sdb, "./Rasters/model-in/sdb-10mmask.tif", format="GTiff", NAflag = NaN)
sdb = brick("./Rasters/model-in/sdb-10mmask.tif")

#Model Points
train.dat = shapefile("./GroundTruth/TrainingSites/ModelInput-09232020-poly/in.dat/FSP-boa-destripe.shp")
train.dat = extract( sdb,train.dat, sp=T)
train.dat = extract( dii,train.dat, sp=T)
train.dat = extract(pca.dat$map,train.dat, sp=T)
shapefile(train.dat, "./GroundTruth/TrainingSites/ModelInput-09232020-poly/in.dat/FSP-boa-destripe-wcc-pca-sdb.shp")
train.dat = shapefile("./GroundTruth/TrainingSites/ModelInput-09232020-poly/in.dat/FSP-boa-destripe-wcc-pca-sdb.shp")
train.dat = train.dat@data

#Add in new VIS points
ab = list.files("./GroundTruth/TrainingSites/ModelInput-09232020-poly/IterAddPoly/",".shp")
for (i in 1:18){
  tmp.b = shapefile(paste0("./GroundTruth/TrainingSites/ModelInput-09232020-poly/IterAddPoly/",ab[i]))
  tmp.b1 = extract( wv.dat,tmp.b, cellnumbers=T,df=T)
  tmp.b1$ndvi = (tmp.b1$re-tmp.b1$r)/(tmp.b1$re+tmp.b1$r)
  tmp.b1$rg = tmp.b1$r/tmp.b1$g
  tmp.b2 = extract( sdb,tmp.b, cellnumbers=T,df=T)
  tmp.b3 = extract( dii,tmp.b, cellnumbers=T,df=T)
  tmp.b4 = extract(pca.dat, tmp.b,cellnumbers=T,df=T)
  tmp.b5 = extract(depth.dat,tmp.b, cellnumbers=T,df=T)
  tmp.b = cbind(tmp.b1[,-c(1)], tmp.b2[,-c(1:2)],tmp.b3[,-c(1:2)],tmp.b4[,-c(1:2)], tmp.b5[,-c(1:2)])
  names(tmp.b)[12] = "sdb_10m"
  names(tmp.b)[29] = "depth"
  tmp.b$InptLbl = 0
  tmp.b$station = strsplit(ab[i], "_")[[1]][1]
  tmp.b$QCC_WV3 = "Pass"
  train.dat = rbind(train.dat[,names(tmp.b)],tmp.b)
 rm(tmp.b,tmp.b1,tmp.b2,tmp.b3,tmp.b4,tmp.b5)}#Bare Poly
for (i in 19:32){
  tmp.b = shapefile(paste0("./GroundTruth/TrainingSites/ModelInput-09232020-poly/IterAddPoly/",ab[i]))
  tmp.b1 = extract( wv.dat,tmp.b, cellnumbers=T,df=T)
  tmp.b1$ndvi = (tmp.b1$re-tmp.b1$r)/(tmp.b1$re+tmp.b1$r)
  tmp.b1$rg = tmp.b1$r/tmp.b1$g
  tmp.b2 = extract( sdb,tmp.b, cellnumbers=T,df=T)
  tmp.b3 = extract( dii,tmp.b, cellnumbers=T,df=T)
  tmp.b4 = extract(pca.dat, tmp.b,cellnumbers=T,df=T)
  tmp.b5 = extract(depth.dat,tmp.b, cellnumbers=T,df=T)
  tmp.b = cbind(tmp.b1[,-c(1)], tmp.b2[,-c(1:2)],tmp.b3[,-c(1:2)],tmp.b4[,-c(1:2)], tmp.b5[,-c(1:2)])
  names(tmp.b)[12] = "sdb_10m"
  names(tmp.b)[29] = "depth"
  tmp.b$InptLbl = 1
  tmp.b$station = strsplit(ab[i], "_")[[1]][1]
  tmp.b$QCC_WV3 = "Pass"
  train.dat = rbind(train.dat[,names(tmp.b)],tmp.b)
  rm(tmp.b,tmp.b1,tmp.b2,tmp.b3,tmp.b4,tmp.b5)}#veg Poly
get.cord = xyFromCell(wv.dat[[1]], train.dat$cell)
out.train.data = SpatialPointsDataFrame(get.cord,train.dat)
out.train.data@proj4string = wv.dat@crs
shapefile(out.train.data, "./GroundTruth/TrainingSites/ModelInput-09232020-poly/in.dat/FSP-VIS-boa-destripe-wcc-pca-sdb.shp")
