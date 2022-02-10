library(dplyr)
library(raster)
library(fields)
library(pracma)
library(ggplot2)
library(ggimage)
library(reshape2)
library(gplots)

##Define all input information here
wv.dat = brick("C://rasterdat.tif")#where is your raster data is to be destriped
names(wv.dat) = c("cb","b","g","y", "r", "re","n1","n2")
dir.out = "C://StripeCorrection//"#directory that must already exist where all output files will be written to
#Define to determine stripe location
#Step 1 values
quant.bpv = 0.85#quantile to define as bright pixels for masking
edge.offset = 113#number of edge columns to exclude from stripe finding, image specific
#Step 2 values
peakheight =  0.0003 #Peak height threshold, image specific
peakdistance = 500 #Peak distance threshold , image specific
lag.n.in = 30 #Lag threshold, image specific
#Step 4 values
lag.n1 = 0:5#how far from the identified stripe start stop to test to flatten the curve for the linear adjustment, image specific
lag.n2 = 20:30#how far from the identified stripe start stop to test to flatten the curve for the linear adjustment, image specific
#Extra values to define
pr = c(0.045,0.0250,0.015,0.020,0.015,0.015,0.005,0.010)#this is just a plotting parameter to make the maps nicer
##
####
 
#Semiautomatic stripe correction once above parameters are defined
#No code needs changing after this line
final.correction.values = matrix(0,ncol=dim(wv.dat)[3],nrow=dim(wv.dat)[2])
colnames(final.correction.values)  =names(wv.dat)
bpv = rep(NA,8)#What is the quantile threshold to mask bright pixels
for(i in 1:dim(wv.dat)[3]){#number of wv bands
  print(i)
xr = as.matrix(wv.dat[[i]])
#Mask out bright pixels
bpv[i] = quantile(na.omit(as.vector(xr)),quant.bpv)#impact of shallow water
xr[xr > bpv[i]] = NaN
# Here I just compute the mean for a column
# and then compute the difference with a lag 
yr = apply(xr,2,"mean",na.rm=T)
diffr = abs(diff(yr[edge.offset:(length(yr)-edge.offset)], lag = lag.n.in ))#Remove edges when doing this
# This is an R function that will find peaks, i.e., jumps in 
# mean reflectance for the image
respeak = findpeaks(diffr, nups = 1, 
                    minpeakdistance = peakdistance ,
                    minpeakheight = peakheight)
if(is.null(respeak)!=T){
respeak[,2:4]=respeak[,2:4]+edge.offset#add edge offset # back in to have correct column index
##
# Here we find the location and magnitude of peaks
peakid = sort(respeak[,2], index.return=T)$ix
peakidx = respeak[peakid,2]
peakamp = respeak[peakid,1]
peak.start = respeak[peakid,3]
peak.stop = respeak[peakid,4]
# The correction of each peak is incremental
coryr = yr
corxr = as.matrix(wv.dat[[i]])
##Which stripe is "widest"
refstripe = c(1,peakidx,length(yr))
refstripe = which.max(diff(refstripe))
#Define correction value per column
coryr.out = yr
coryr.out = ifelse(coryr.out>0,0,coryr.out)
#Generate correction values "right" of the stripe
cor.value.left = rep(NA, length(peakidx[c(refstripe:(length(peakidx)))]))
use.index = c(1,peakidx,length(yr)+1)
for (j in 1:length(cor.value.left)){
  dim.use = j+refstripe
  #don't use absolute because want to decrease if it is negative
  cor.value.left[j] =   (mean(coryr[(use.index[dim.use]-52):(use.index[dim.use] - 42)]) -
                           mean(coryr[(use.index[dim.use]+42):(use.index[dim.use] + 52)]) )
  coryr[use.index[dim.use]:(use.index[dim.use+1]-1)] = coryr[use.index[dim.use]:(use.index[dim.use+1]-1)] + cor.value.left[j]
  coryr.out[use.index[dim.use]:(use.index[dim.use+1]-1)] = coryr.out[use.index[dim.use]:(use.index[dim.use+1]-1)] + cor.value.left[j]
  corxr[,use.index[dim.use]:(use.index[dim.use+1]-1)] = corxr[,use.index[dim.use]:(use.index[dim.use+1]-1)] + cor.value.left[j]}
#Correct "left" of the stripe
if(refstripe>1){
dim.use = seq(refstripe,2,-1)
cor.value.right = rep(NA,length(dim.use))
for(j in 1:length(dim.use)){
  cor.value.right[j] =   (mean(coryr[(use.index[dim.use[j]]+52):(use.index[dim.use[j]] + 42)]) -
                            mean(coryr[(use.index[dim.use[j]]-42):(use.index[dim.use[j]] - 52)]) )
  coryr[use.index[dim.use[j]-1]:(use.index[dim.use[j]]-1)] = coryr[use.index[dim.use[j]-1]:(use.index[dim.use[j]]-1)] + cor.value.right[j]
  coryr.out[use.index[dim.use[j]-1]:(use.index[dim.use[j]]-1)] = coryr.out[use.index[dim.use[j]-1]:(use.index[dim.use[j]]-1)] + cor.value.right[j]
  corxr[,use.index[dim.use[j]-1]:use.index[dim.use[j]]] = corxr[,use.index[dim.use[j]-1]:use.index[dim.use[j]]] + cor.value.right[j]}}
#

##Output column correction data
csv.coryr = cbind(yr,coryr,coryr.out)
if(min(csv.coryr[,2],na.rm = T)<0){csv.coryr[,3] = csv.coryr[,1] - csv.coryr[,2]}
colnames(csv.coryr) = c("ColMean","ColMeanCor","ColCorValue")
write.csv(csv.coryr,paste0(dir.out,i,"-bpv",quant.bpv,"-ph",peakheight ,"-pd",peakdistance,"-ln",lag.n.in , "-corval.csv"),row.names = F)
### Making PNG...
png(paste0(dir.out,i,"-bpv",quant.bpv,"-ph",peakheight ,"-pd",peakdistance,"-ln",lag.n.in , ".png"),width=6.5,height=6,units="in",res=300)
par(mfcol=c(3,1),mar=c(2,3,0,0),oma=c(0,0,0.5,0.5))
pr2=c(seq(minValue(wv.dat[[i]]),pr[i],(pr[i]-minValue(wv.dat[[i]]))/10),maxValue(wv.dat[[i]]))
plot(raster(as.matrix(wv.dat[[i]])),maxpixels=(dim(xr)[1]*dim(xr)[1])/100,legend=T,
     col=rich.colors(length(pr2)-1),useRaster=F, breaks=pr2 )
text(x=0.1,y=0.9, paste0("band-", i,"-uncorrected"))
plot(raster(corxr),maxpixels=(dim(xr)[1]*dim(xr)[1])/100,legend=T,
     col=rich.colors(length(pr2)-1),useRaster=F, breaks=pr2)
text(x=0.1,y=0.9, paste0("band-", i,"-corrected"))
plot(coryr,type="l",xaxs="i",col="red",ylim=c(pr[i]-0.015,pr[i]) )
lines (yr,type="l",xaxs="i")
abline(v=respeak[,4],col=3)#,xpd=NA
legend ("topleft", title=i,legend = c("Uncorrected","Corrected",
                   paste0("bpv=",quant.bpv), paste0("peakheight=",peakheight), paste0("peakdistance=",peakdistance), 
                   paste0("lag=",lag.n.in )), col=c("black","red",rep("NA",4)),lty=1, bty="n",ncol=2)
dev.off()
if(exists("cor.value.right")==T){
  data.out = cbind(peakamp, peakidx,peak.start,peak.stop,rep(refstripe,length(peakamp)),c(rev(cor.value.right),cor.value.left))}
if(exists("cor.value.right")==F){
  data.out = cbind(peakamp, peakidx,peak.start,peak.stop,rep(refstripe,length(peakamp)),c(cor.value.left))}
#output respeak data
colnames(data.out)= c("amp","id","start","stop","refstripe","newcorrectionvalue")
write.csv(data.out,paste0(dir.out,i,"-bpv",quant.bpv,"-ph",peakheight ,"-pd",peakdistance,"-ln",lag.n.in , ".csv"),row.names = F)
rm(corxr,xr,cor.value.left,coryr,diffr,dim.use,j,peakamp,peakid,peakidx,refstripe,use.index,yr,
   peak.start,peak.stop,data.out,csv.coryr,coryr.out)
if(exists("cor.value.right")==T){rm(cor.value.right)}}
if(is.null(respeak)==T){
  write.csv(1,paste0(dir.out,i,"-bpv",quant.bpv,"-ph",peakheight ,"-pd",peakdistance,"-ln",lag.n.in , "-corval.csv"),row.names = F)
  write.csv(1,paste0(dir.out,i,"-bpv",quant.bpv,"-ph",peakheight ,"-pd",peakdistance,"-ln",lag.n.in , ".csv"),row.names = F)
  ### Making PNG...
  png(paste0(dir.out,i,"-bpv",quant.bpv,"-ph",peakheight ,"-pd",peakdistance,"-ln",lag.n.in , ".png"),width=6.5,height=6,units="in",res=300)
  par(mfcol=c(3,1),mar=c(2,3,0,0),oma=c(0,0,0.5,0.5))
  pr2=c(seq(minValue(wv.dat[[i]]),pr[i],(pr[i]-minValue(wv.dat[[i]]))/10),maxValue(wv.dat[[i]]))
  plot(raster(as.matrix(wv.dat[[i]])),maxpixels=(dim(wv.dat)[1]*dim(wv.dat)[1])/100,legend=T,
       col=rich.colors(length(pr2)-1),useRaster=F, breaks=pr2 )
  text(x=0.1,y=0.9, paste0("band-", i,"-uncorrected"))
  plot(raster(as.matrix(wv.dat[[i]])),maxpixels=(dim(wv.dat)[1]*dim(wv.dat)[1])/100,legend=T,
       col=rich.colors(length(pr2)-1),useRaster=F, breaks=pr2)
  text(x=0.1,y=0.9, paste0("band-", i,"-corrected"))
  plot(yr,type="l",xaxs="i",col="red",ylim=c(pr[i]-0.015,pr[i]) )
   dev.off()
  
}
rm(respeak)
}
#Jump correction, step 4 in manuscript
cor.files = list.files(dir.out, "corval.csv")#output in loop 1
respeak.files = list.files(dir.out, "30.csv")#output in loop 1
for(i in 1:dim(wv.dat)[3]){
  cor.values = read.csv(paste0(dir.out,cor.files[i]))
  respeak = read.csv(paste0(dir.out,respeak.files[i]))
  if(dim(respeak)[2]>1){
  lag.n = expand.grid(lag.n1,lag.n2)
  lag.n.min = rep(NA, length(lag.n[,1]))
  out.dim.lag = rep(NA, length(respeak$start))
  
  for(j in 1:length(respeak$start)){
    for(k in 1:length(lag.n[,1])){
      coryr.out.ln = cor.values$ColCorValue
      #Get start and stop x values including extras values past if start/stop was not sensitive enough
      x.dat = c(respeak$start[j]-lag.n[k,1],respeak$stop[j]+lag.n[k,2])
      #Get start and stop correction values
      y.dat = c(coryr.out.ln[respeak$start[j]-lag.n[k,1]],coryr.out.ln[respeak$stop[j]+lag.n[k,2]])
      #Linearlly interpolate new correction values between start and stop
      lm.out = lm(y.dat~x.dat)
      new.x.dat = seq(x.dat[1]+1,x.dat[2]-1,1)
      new.y.dat = (new.x.dat *lm.out$coefficients[2])+lm.out$coefficients[1]
      coryr.out.ln[new.x.dat] = new.y.dat
      #find the difference between the max and minimum values to make transition as flat as possible
      #avoid slope as if there are two opposite peaks slope is zero still
      lag.n.min[k]= max((cor.values$ColMean+coryr.out.ln)[new.x.dat])-min((cor.values$ColMean+coryr.out.ln)[new.x.dat])
    }
    #Write out the dimension of lag.n to use per stripe
    out.dim.lag[j] =which.min(lag.n.min)
  }
  ##Identify and correct  the "peaks" that occur over the peaks
  png(paste0(dir.out,"band",i,".png"),width=7.5,height=7.5,units="in",res=300)
  par(mfcol=c(ceiling(sqrt(length(respeak$start))),ceiling(sqrt(length(respeak$start)))),mar=c(3,3,1,1))
  coryr.out.ln = cor.values$ColCorValue
  coryr.out.ln.buffer = cor.values$ColCorValue
  for(j in 1:length(respeak$start)){
    plot(cor.values$ColMean ,type="b",xaxs="i",col="black",pch=20,xlim =c(respeak$start[j]-75,respeak$stop[j]+75),
         ylim=c(min(cor.values$ColMean[respeak$start[j]-75:respeak$stop[j]+75])-0.0015,
               max(cor.values$ColMean[respeak$start[j]-75:respeak$stop[j]+75])+0.0015))
    lines(cor.values$ColMeanCor,col="red",type="b",pch=1)
    #correct without adding buffer points
    x.dat = c(respeak$start[j],respeak$stop[j])
    y.dat = c(coryr.out.ln[respeak$start[j]],coryr.out.ln[respeak$stop[j]])
    lm.out = lm(y.dat~x.dat)
    new.x.dat = seq(x.dat[1]+1,x.dat[2]-1,1)
    new.y.dat = (new.x.dat *lm.out$coefficients[2])+lm.out$coefficients[1]
    coryr.out.ln[new.x.dat] = new.y.dat
    if ( min(cor.values$ColMean,na.rm = T)<0) {lines((cor.values$ColMean-coryr.out.ln),col="blue",type="b",pch=20)}
    if ( min(cor.values$ColMean,na.rm = T)>0) {lines((cor.values$ColMean+coryr.out.ln),col="blue",type="b",pch=20)}
    #correct adding buffer points
    x.dat = c(respeak$start[j]-lag.n[out.dim.lag[j],1],respeak$stop[j]+lag.n[out.dim.lag[j],2])
    y.dat = c(coryr.out.ln.buffer[respeak$start[j]-lag.n[out.dim.lag[j],1]],coryr.out.ln.buffer[respeak$stop[j]+lag.n[out.dim.lag[j],2]])
    lm.out = lm(y.dat~x.dat)
    new.x.dat = seq(x.dat[1]+1,x.dat[2]-1,1)
    new.y.dat = (new.x.dat *lm.out$coefficients[2])+lm.out$coefficients[1]
    coryr.out.ln.buffer[new.x.dat] = new.y.dat
    if ( min(cor.values$ColMean,na.rm = T)<0){lines((cor.values$ColMean-coryr.out.ln.buffer),col="green",type="b",pch=20)}
    if ( min(cor.values$ColMean,na.rm = T)>0){lines((cor.values$ColMean+coryr.out.ln.buffer),col="green",type="b",pch=20)}
    #show where correction was extended to
    abline(v= (respeak$start[j]-lag.n[out.dim.lag[j],1]))#,xpd=NA)
    abline(v=(respeak$stop[j]+lag.n[out.dim.lag[j],2]))#,xpd=NA)
    legend("bottomright",pch=c(20,1,20,20),col=c("black","red","blue","green"),lty=1,bty="n",ncol=2,
           title = ifelse(out.dim.lag[j]>1,"","B=G" ),
           legend=c("Uncor","Cor","Lm Cor","Lm Cor+"))
  }
  dev.off()
  final.correction.values[,i] =  coryr.out.ln.buffer
}}
write.csv(final.correction.values, paste0(dir.out, "allband-bestcorval-allband.csv"),row.names = F)
#Write Raster Out
for(i in 1:8){
  cor.matrix.dat = matrix(final.correction.values[,i],nrow=dim(wv.dat[1])) 
  cor.matrix = raster(cor.matrix.dat)
  crs(cor.matrix) =crs(wv.dat)
  extent(cor.matrix) =extent(wv.dat)
  res(cor.matrix) = res(wv.dat)
  cor.matrix[] = cor.matrix.dat
  ##
  cor.values = read.csv(paste0(dir.out,cor.files[i]))
  respeak = read.csv(paste0(dir.out,respeak.files[i]))
  ##Plotting
  ### Making PNG...
  png(paste0(dir.out,i,"-bpv",quant.bpv,"-ph",peakheight ,"-pd",peakdistance,"-ln",
             lag.n.in , "-jumpcrc.png"),width=6.5,height=6,units="in",res=300)
  par(mfcol=c(3,1),mar=c(2,3,0,0),oma=c(0,0,0.5,0.5))
  pr2=c(seq(minValue(wv.dat[[i]]),pr[i],(pr[i]-minValue(wv.dat[[i]]))/10),maxValue(wv.dat[[i]]))
  plot(raster(as.matrix(wv.dat[[i]])),maxpixels=(dim(wv.dat)[1]*dim(wv.dat)[1])/100,legend=T,
       col=rich.colors(length(pr2)-1),useRaster=F, breaks=pr2 )
  text(x=0.1,y=0.9, paste0("band-", i,"-uncorrected"))
  ##
  if ( min(cor.values$ColMean,na.rm = T)<0) {wv.dat[[i]] = wv.dat[[i]]-cor.matrix}
  if ( min(cor.values$ColMean,na.rm = T)>0) {wv.dat[[i]] = wv.dat[[i]]+cor.matrix}
  plot(raster(as.matrix(wv.dat[[i]])),maxpixels=(dim(wv.dat)[1]*dim(wv.dat)[1])/100,legend=T,
       col=rich.colors(length(pr2)-1),useRaster=F, breaks=pr2)
  text(x=0.1,y=0.9, paste0("band-", i,"-Finalcorrected"))
  ##
 
  if(dim(respeak)[2]>1){
  plot(cor.values$ColMean ,type="l",xaxs="i",col="black",ylim=c(pr[i]-0.015,pr[i])  )
  lines (cor.values$ColMeanCor,type="l",xaxs="i",col="red")
  if ( min(cor.values$ColMean,na.rm = T)<0) {lines(cor.values$ColMean-final.correction.values[,i],col="green",type="l",xaxs="i",lty=2)}
  if ( min(cor.values$ColMean,na.rm = T)>0) {lines(cor.values$ColMean+final.correction.values[,i],col="green",type="l",xaxs="i",lty=2)}
  abline(v=respeak[,4],col="grey")#,xpd=NA
  legend ("topleft", title=i,legend = c("Uncorrected","Corrected","CorrectedFinal"), col=c("black","red","green"),
          lty=c(1,1,2), bty="n",ncol=3)}
  dev.off()
  ##
  
  rm(cor.matrix.dat,cor.matrix)
}
writeRaster(wv.dat,paste0(dir.out, "BOA-destripe.tif"),format="GTiff",NAflag = NaN,overwrite=T)#this is your final destriped data
##