rm(list = ls())
library("raster")
library("rgdal")

setwd("C://Users//wilsonkri//Documents//Backup-WV//20190817//")

##
depth= raster("./Extrafiles/Depth.tif")
raster.dat = "./StripeCorrection/PD3-use//BOA-destripe.tif"
raster.name.dat = c("cb","b","g","y","r","re","n1","n2")
raster.dim.use = 5 #number of bands to use the WCC
depth.dat = "./Extrafiles/Depth.tif"
write.data = "./WCC_Destripe/Lyzenga1985-ds-alldpth.tif"
write.folder = "./WCC_Destripe/"
shape.dat = "./WCC_BOA/Sand/"

#Read data
dat = stack(raster.dat,depth.dat)
names(dat) = c(raster.name.dat,"depth")

#Subset polygon data
files = list.files(shape.dat, pattern=".shp")
shape.in = shapefile(paste(shape.dat,files[1],sep=""))
shape.in@data$pol = 1
for (i in 2:length(files)){
  tmp = shapefile(paste(shape.dat,files[i],sep=""))
  tmp@data$pol = i
  shape.in =bind(shape.in,tmp)
  rm(tmp)
}
rm(files,i)
shape.in = spTransform(shape.in, dat@crs)
plot(dat[[1]])
plot(shape.in,add=T)
dat = extract(x=dat,y=shape.in,df=T)
dat = dat[,c("ID",raster.name.dat[1:raster.dim.use],"depth")]
rm(shape.in)

#Generate ln
dat =na.omit(dat)
dat.ln=dat
dat.ln[,c(raster.name.dat[1:raster.dim.use])] = log (dat[,c(raster.name.dat[1:raster.dim.use])])#take natural logarithm
#Calculate variance
sub.dat.variance = apply(dat.ln, 2, var) #calculate variance 
#Generate unique combinations
val.in = expand.grid(names(sub.dat.variance[2]),names(sub.dat.variance[3:(1+raster.dim.use)]))
i=3
{while(i<=raster.dim.use){
    print(i)
    test = expand.grid(names(sub.dat.variance[i]),names(sub.dat.variance[(i+1):(1+raster.dim.use)]))
    val.in = rbind(val.in,test)
  rm(test)
  i = i+1}}#generate unqiue DDI
val.in = cbind(val.in,rep(NA,length(val.in[,1])),rep(NA,length(val.in[,1])),rep(NA,length(val.in[,1]))
               ,rep(NA,length(val.in[,1])),rep(NA,length(val.in[,1])),rep(NA,length(val.in[,1])))
names(val.in) = c("B1","B2","Variance.B1","Variance.B2", "Covariance","a","slope","R2")
val.in[,1] = as.character(val.in[,1])
val.in[,2] = as.character(val.in[,2])
rm(i)
#Plot Original Values
for ( i in 1:length(val.in[,1])){
  png(paste(write.folder,"biplot-",val.in[i,1],val.in[i,2],"-ori.png",sep=""),pointsize=10,family="serif")
  plot(dat[,val.in[i,1]], dat[,val.in[i,2]], ylab = paste (val.in[i,2],sep=""),
       xlab = paste (val.in[i,1],sep=""),pch=19, col=dat.ln[,1])
  test = lm(dat[,val.in[i,2]]~ dat[,val.in[i,1]]) 
  abline(a=  test$coefficients[1],b=test$coefficients[2])
  legend("bottomright", legend=paste0("R2=", round(summary(test)$r.squared,3)), bty="n")
  dev.off()
   rm(test)}
#Calculate covariance
for ( i in 1:length(val.in[,1])){
  val.in$Covariance[i] = cov(dat.ln[,val.in[i,1]], dat.ln[,val.in[i,2]])
}
#Put variance into val.in
i=1
{while(i<raster.dim.use){
  i = i+1
  test = ifelse(names(sub.dat.variance[i])==val.in$B1,sub.dat.variance[i],NA)
  val.in$Variance.B1 = ifelse(is.na(test)==T,val.in$Variance.B1,test)
}}
i=1
{while(i<=raster.dim.use){
  i = i+1
  test = ifelse(names(sub.dat.variance[i])==val.in$B2,sub.dat.variance[i],NA)
  val.in$Variance.B2 = ifelse(is.na(test)==T,val.in$Variance.B2,test)
  rm(test)
}}
#calculate a
for ( i in 1:length(val.in[,1])){
  val.in$a[i] = (val.in$Variance.B1[i]-val.in$Variance.B2[i])/(2*val.in$Covariance[i]) 
}
#calculate the slope
for ( i in 1:length(val.in[,1])){
  val.in$slope[i] = val.in$a[i] + sqrt((val.in$a[i]^2)+1)
}
#Plot each bi-plot
for ( i in 1:length(val.in[,1])){
  png(paste(write.folder,"biplot-",val.in[i,1],val.in[i,2],".png",sep=""),pointsize=10,family="serif")
  plot(dat.ln[,val.in[i,1]], dat.ln[,val.in[i,2]], ylab = paste ("ln ", val.in[i,2],sep=""),
       xlab = paste ("ln ", val.in[i,1],sep=""),pch=19, col=dat.ln[,1])
  test = lm(dat.ln[,val.in[i,2]]~ dat.ln[,val.in[i,1]]) 
  abline(a=  test$coefficients[1],b=test$coefficients[2])
  legend("bottomright", legend=paste0("R2=", round(summary(test)$r.squared,3)), bty="n")
  dev.off()
  val.in$R2[i] = round(summary(test)$r.squared,3)
  print(val.in[i,])
  print(test)
  rm(test)}
#Write out the coefficients
write.csv(val.in,paste(write.folder,"WCC-Coef",".csv",sep="") )
##
rm(dat,dat.ln,depth.dat,sub.dat.variance,i)

#Apply the  water column correction to the visible bands
ras.dat = brick(raster.dat)
names(ras.dat) = raster.name.dat
ras.dat = ras.dat[[1:raster.dim.use]]
ras.dat = log(ras.dat)
names(ras.dat) = raster.name.dat[1:raster.dim.use]
dii = stack(ras.dat,ras.dat)
names(dii) = c(paste0("dii",val.in[,1],val.in[,2]))
dii=dii*0
for (i in 1:length(val.in[,1])){
  print(i)
  dii[[i]] = ras.dat[[val.in[i,1]]]-(val.in$slope[i]*ras.dat[[val.in[i,2]]])
}
#Export data
rm(ras.dat)
names(dii) = c(paste0("dii",val.in[,1],val.in[,2]))
writeRaster(dii, write.data, format="GTiff", NAflag = NaN)
###




