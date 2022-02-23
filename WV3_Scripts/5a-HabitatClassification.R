library("raster")
library("caret")
library("irr")
library("doParallel")  #Foreach Parallel Adaptor 
library("foreach")     #Provides foreach looping construct
rasterOptions(maxmemory=7e+09,chunksize = 4e+08)

setwd("C:\\Where\\is\\my\\data\\")

#Read in raster dataset
wv.dat = stack("BOA-destripe.tif",
               "Lyzenga1985-ds-alldpth.tif",
               "PCA-nostand-Destripe-Level2-reflectance-20190817.tif")
names(wv.dat) = c( "cb", "b", "g" , "y"  , "r" ,"re" , "n1", "n2" ,        
                   "dii_cbb","dii_cbg","dii_cby","dii_cbr","dii_bg","dii_by", "dii_br" ,"dii_gy", "dii_gr"  , "dii_yr",           
                    "PCA1","PCA2",  "PCA3", "PCA4", "PCA5" , "PCA6" )
depth.mask=wv.dat[[1]]#only used for blank raster to fill
#dep.val = 10 
ndvi.thres = 0.3 
b.thres = 0.03 
gr = 0.42
#Pick one training dataset
train.dat = shapefile("AllPolyDS.shp")
names(train.dat)[4] = "InptLbl"
#
list.bands.in = list( c("InptLbl",     "PCA1","PCA2",  "PCA3","PCA4"),
                      c("InptLbl",     "PCA1","PCA2",  "PCA3","PCA4","PCA5"),
                      c("InptLbl",     "PCA1","PCA2",  "PCA3"),
                      c("InptLbl",     "PCA1","PCA2"),
                      c("InptLbl",  "cb","b","g","y","r","re"),
                      c("InptLbl",  "b","g","r"),
                      c("InptLbl",  "cb","b","g","y","r"),
                      c("InptLbl",  "cb","b","g","y"),
                      c("InptLbl", "b","g","y","r","re"),
                      c("InptLbl", "b","g","y","r"),
                      c("InptLbl", "b","g","y"),
                      c("InptLbl","dii_cbb","dii_cbg","dii_cby","dii_cbr","dii_bg","dii_by","dii_br","dii_gy","dii_gr","dii_yr"),
                      c("InptLbl",     "PCA1","PCA2",  "PCA3","PCA4",   
                             "dii_cbb","dii_cbg","dii_cby","dii_cbr","dii_bg","dii_by","dii_br","dii_gy","dii_gr","dii_yr",  
                             "cb","b","g","y","r","re"))
list.out.folder = c("./RF/")
UseCores = detectCores() -1-5-5-1


##Generate Training Data 
train.dat2 = train.dat
train.dat = train.dat@data
train.dat$InptLbl = as.numeric(train.dat$InptLbl)
train.dat$InptLbl = ifelse(train.dat$InptLbl==0,0,
                           ifelse(train.dat$InptLbl<0,NA,#Omit if keeping ODW
                           ifelse(is.na(train.dat$InptLbl)==T,NA,1)))
train.dat = train.dat[!is.na(train.dat$InptLbl),]
train.dat = train.dat[!is.na(train.dat$b),]
train.dat$ndvi = (train.dat$re - train.dat$r)/(train.dat$re + train.dat$r)
train.dat$rg = train.dat$r/train.dat$g
train.dat = train.dat[train.dat$Depth<7, ]
##
wv.dat=getValues(wv.dat[[c( "cb", "b", "g" , "y"  , "r" ,"re" ,         
                            "dii_cbb","dii_cbg","dii_cby","dii_cbr","dii_bg","dii_by", "dii_br" ,"dii_gy", "dii_gr"  , "dii_yr",           
                            "PCA1","PCA2",  "PCA3", "PCA4", "PCA5" , "PCA6" )]])
id.na = which(!is.na(wv.dat[,1]))#gives index of true values
wv.dat = na.omit(wv.dat)#remove NA
wv.dat=as.data.frame(wv.dat)
##
wv.dat$ndvi = (wv.dat$re-wv.dat$r)/(wv.dat$re+wv.dat$r)
wv.dat$rg = wv.dat$r/wv.dat$g
gc()
#

for (k in 1:length(list.bands.in)){
#
bands.in = list.bands.in[k]
bands.in = unlist(bands.in)
out.folder = paste(bands.in, collapse = "-")
out.folder = strsplit(out.folder, "InptLbl-" )[[1]][2]
out.folder = paste0(list.out.folder,out.folder,"/")
#

##Create directoriers
dir.create(out.folder)
dir.create(paste(out.folder,"/cvrun",sep=""))
dir.create(paste(out.folder,"/cvrunmf",sep=""))

#Define Model tuning
ras.value = train.dat[,bands.in]
ras.value[,1] = as.factor(ras.value[,1])#dependent value as a factor
rfFit = train(form = InptLbl ~ . , data = ras.value,
              method = "rf", tuneLength = (dim(ras.value)[2]-2),
              trControl = trainControl(method = "repeatedcv",number = 5,repeats = 10),
              verbose = TRUE)
saveRDS(rfFit, paste(out.folder,"tunemodel",sep=""))
##
##Define model with all training data and best tune
num.mtry = which.max(rfFit$results$Accuracy)+1
num.mtry = expand.grid(mtry = num.mtry)#only use best mtry evaulated above
rm(rfFit)
#
set.seed(805)
trainIndex = createMultiFolds(ras.value[,1], k = 5, times = 10)#generate mutliplte data splits
###
registerDoSEQ()#removes previous clusters
#Define how many cores you want to use and Register CoreCluster
cl  = makeCluster(UseCores)
registerDoParallel(cl) 
#
foreach(i=1:50) %dopar% {
  library("raster")
  library("caret")
  library("irr")
  #Out data
  full.model = paste(out.folder, "cvrun/", "k",i, ".rds",sep="")
  exp.file = paste(out.folder, "cvrun/", "k",i,".csv",sep="")
  out.ras = paste(out.folder, "cvrun/","k",i,".tif",sep="")
  out.ras1 = paste(out.folder, "cvrunmf/","k",i,"-mf.tif",sep="")
  exp.file1 = paste(out.folder, "cvrunmf/", "k",i,".csv",sep="")
  #Build model
  ras.valueTrain = ras.value[ trainIndex[[i]],]
  ras.valueTest  = ras.value[-trainIndex[[i]],]
  rfFit.k = train(form = InptLbl ~ . , data = ras.valueTrain,
                  method = "rf", 
                  tuneGrid = data.frame(num.mtry),
                  trControl = trainControl(method = "none"),
                  verbose = FALSE)
  #
  saveRDS(rfFit.k$finalModel, full.model)
  ##
  set.seed(3)
  rfFit.test = predict(rfFit.k, ras.valueTest)
  #fix threshold based points
  rfFit.threshold = train.dat[-trainIndex[[i]],]
  rfFit.final = ifelse(rfFit.threshold$ndvi>=ndvi.thres , 2,
                       ifelse(rfFit.threshold$ndvi<ndvi.thres & rfFit.threshold$b >= b.thres,1,
                              ifelse(rfFit.threshold$ndvi<ndvi.thres & rfFit.threshold$rg <= gr, 1,
                                       rfFit.test)))
  #Evaulate model
  full.dat=as.data.frame(cbind(rfFit.final,ras.valueTest[,1]))
  conmat = table(full.dat)
  conmat = apply(conmat, 2, function(x) as.numeric(as.character(x)))
  conkapp = kappa2(full.dat, weight = c("unweighted"))
  n = sum(conmat)#number of observations
  oa = sum(diag(conmat))#number of correct classifications overall 
  oa.per = oa/n*100#percentage of correct classifications overall
  colsums = apply(conmat, 2, sum)
  PA = diag(conmat) / colsums*100
  rowsums = apply(conmat, 1, sum)
  UA = diag(conmat) / rowsums*100
  export = rbind(conmat,PA)
  export = as.matrix(cbind(export,c(UA,oa.per),c(conkapp$value,conkapp$p.value,rep(NA,nrow(export)-2))))
  write.csv(export,exp.file, row.names = F)
  rm(export, full.dat,conmat,conkapp,n,oa,oa.per,colsums,PA,rowsums,UA)
  ##
  set.seed(3)
  out.dat = predict(rfFit.k, wv.dat)
  ##
  out.dat = as.numeric(out.dat)
  out.dat = out.dat-1
  out.dat= ifelse(wv.dat[,"ndvi"]>=ndvi.thres, 1,
                  ifelse(wv.dat[,"ndvi"]<ndvi.thres & wv.dat[,"b"] >= b.thres,0,
                         ifelse(wv.dat[,"ndvi"]<ndvi.thres & wv.dat[,"rg"] <= gr, 0,
                                 out.dat)))
  ##
  out.dat.raster = raster(depth.mask)#define a raster
  out.dat = as.integer(out.dat)#define as integer to reduce file size
  out.dat.raster[id.na] =  out.dat#save output
  writeRaster(out.dat.raster, out.ras, format="GTiff",NAflag = NaN,overwrite=T)
  ##Add focal filter
  mj.raster = focal(out.dat.raster, w=matrix(1,3,3), modal,na.rm=T)
  #mj.raster = raster(out.ras1)
  mj.raster = mask(x=mj.raster, mask=depth.mask)
  writeRaster(mj.raster, out.ras1, format="GTiff",NAflag = NaN,overwrite=T)
  mj.test = extract(mj.raster,  train.dat[-trainIndex[[i]],"cell"])
  full.dat=as.data.frame(cbind( mj.test,ras.valueTest[,1]))
  full.dat$V2 = ifelse(full.dat$V2==1,0,1)
  conmat = table(full.dat)
  conmat = apply(conmat, 2, function(x) as.numeric(as.character(x)))
  conkapp = kappa2(full.dat, weight = c("unweighted"))
  n = sum(conmat)#number of observations
  oa = sum(diag(conmat))#number of correct classifications overall 
  oa.per = oa/n*100#percentage of correct classifications overall
  colsums = apply(conmat, 2, sum)
  PA = diag(conmat) / colsums*100
  rowsums = apply(conmat, 1, sum)
  UA = diag(conmat) / rowsums*100
  export = rbind(conmat,PA)#
  export = as.matrix(cbind(export,c(UA,oa.per),c(conkapp$value,conkapp$p.value,rep(NA,nrow(export)-2))))
  write.csv(export,exp.file1, row.names = F)
  rm(export, full.dat,conmat,conkapp,n,oa,oa.per,colsums,PA,rowsums,UA)
  ##
  rm(ras.valueTest,ras.valueTrain,rfFit.test,exp.file,out.dat, out.ras, out.dat.raster)
    rm(rfFit.k, full.model)
  }
#end cluster
stopCluster(cl)
registerDoSEQ()#removes previous clusters
rm(cl)
gc()
###
rm(num.mtry,trainIndex,ras.value)

#Generate Mean confusion matrix statistics
a = list.files(path=paste(out.folder, "cvrun/", sep=""), pattern=".csv")
ab = paste(out.folder, "cvrun/",a,sep="")
my.list = lapply(ab, read.csv)
my.list  = lapply(my.list , as.matrix)
rm(a,ab)
mat_rows =nrow(my.list[[1]])
len_mat = length(my.list[[1]])
vec = sapply(1:len_mat, function(j) mean(sapply(1:length(my.list), function(i) as.numeric(as.matrix(my.list[[i]])[j])), na.rm=TRUE))
final_mat = matrix(as.numeric(vec), nrow=mat_rows)
write.csv(final_mat,paste(out.folder, "CMsummary.csv", sep=""))
rm(mat_rows,len_mat,vec,final_mat,my.list)
##

#Generate Mean confusion matrix statistics2
a = list.files(path=paste(out.folder, "cvrunmf/", sep=""), pattern=".csv")
ab = paste(out.folder, "cvrunmf/",a,sep="")
my.list = lapply(ab, read.csv)
my.list  = lapply(my.list , as.matrix)
rm(a,ab)
mat_rows =nrow(my.list[[1]])
len_mat = length(my.list[[1]])
vec = sapply(1:len_mat, function(j) mean(sapply(1:length(my.list), function(i) as.numeric(as.matrix(my.list[[i]])[j])), na.rm=TRUE))
final_mat = matrix(as.numeric(vec), nrow=mat_rows)
write.csv(final_mat,paste(out.folder, "CMsummary-mf.csv", sep=""))
rm(mat_rows,len_mat,vec,final_mat,my.list)
##

##Read in all raster
a = list.files(path=paste(out.folder, "cvrun/", sep=""), pattern=".tif$")
ab = paste(out.folder, "cvrun/",a,sep="")
clas.dat = stack(ab)
rm(a,ab)
#
beginCluster()
#Probability
a = dim(clas.dat)[3]
get.prob = function(x) {round((sum(x)/a)*100,2)}
prob.class.dat = clusterR(clas.dat, calc, args=list(fun=get.prob),export='a')#Generate vegetated probability
writeRaster(prob.class.dat,paste(out.folder, "probability.tif", sep=""),
            format="GTiff",NAflag = NaN,overwrite=T) 
endCluster()
rm(a,clas.dat)

##Read in all raster2
a = list.files(path=paste(out.folder, "cvrunmf/", sep=""), pattern=".tif$")
ab = paste(out.folder, "cvrunmf/",a,sep="")
clas.dat = stack(ab)
rm(a,ab)
#
beginCluster()
#Probability
a = dim(clas.dat)[3]
get.prob = function(x) {round((sum(x)/a)*100,2)}
prob.class.dat = clusterR(clas.dat, calc, args=list(fun=get.prob),export='a')#Generate vegetated probability
writeRaster(prob.class.dat,paste(out.folder, "probabilitymf.tif", sep=""),
            format="GTiff",NAflag = NaN,overwrite=T) 
endCluster()
rm(a,clas.dat)

##
class.dat = extract(prob.class.dat,train.dat2,df=T)
oa.all = rep(NA, 100)
kapp.all = rep(NA, 100)
for (i in 1:100){
  rf.thresh.lab = ifelse((class.dat[,2])>=i,1,0)
  full.dat = as.data.frame(cbind(rf.thresh.lab,train.dat$InptLbl ))
  #Evaulate model
  conmat = table(full.dat)
  conkapp = kappa2(full.dat, weight = c("unweighted"))
  n = sum(conmat)#number of observations
  oa = sum(diag(conmat))#number of correct classifications overall 
  oa.per = oa/n*100#percentage of correct classifications overall
  colsums = apply(conmat, 2, sum)
  PA = diag(conmat) / colsums*100
  rowsums = apply(conmat, 1, sum)
  UA = diag(conmat) / rowsums*100
  export = rbind(conmat,PA)
  export = as.data.frame(cbind(export,c(UA,oa.per),c(conkapp$value,conkapp$p.value,rep(NA,nrow(export)-2))))
  names(export) = c("bare","vegetated","UA","kappa-s-p")
  row.names(export) = c("bare","vegetated","PA")
  oa.all[i]=oa.per
  kapp.all[i] = conkapp$value
  rm(full.dat,conmat,conkapp,n,oa,oa.per,colsums,PA,rowsums,UA,rf.thresh.lab)
}
#write.csv(export,exp.file, row.names = T)
png(paste(out.folder, "threshold.png", sep=""))
par(mfrow=c(1,2),mai=c(0.8,0.8,0.1,0.1))
plot(1:100,oa.all,ylim=c(70,100), xlab=c("Veg Threshold"), ylab=c("OVerall Map Accuracy"),pch=20)
plot(1:100,kapp.all,ylim=c(0.4,1), xlab=c("Veg Threshold"), ylab=c("Kappa"),pch=20)
dev.off()
###
rm(list = ls()[!ls() %in% c("wv.dat", "dep.val", "ndvi.thres","b.thres","gr","train.dat","train.dat2","depth.mask","id.na",
                            "list.bands.in","list.out.folder","UseCores")])
removeTmpFiles(h=0)
gc()

}
