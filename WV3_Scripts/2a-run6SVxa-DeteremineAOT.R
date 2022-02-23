#Use this code to define the best AOT to use
#Test on different bands and polygons
#this requires 6SV to be installed on your computer

##
library(readxl)
#
#Point to the folder where 6SV is installed our your computer
setwd("D:\\SixSV")#For R
shell("cd D:\\SixSV")#For command line script
#
#Pre extract TOA data for various bands
#this is a vector of optically deep water over which you define the AOT
nir.ref = as.vector(read.table("./OpticallyDeepWater.txt",header=T))#this is text file for one band for ODW to define the AOT for
#what is the band number you are using
band.num=8
#WV-3 filter functions must be acquired from the data provider
dir.relative.response = "./WV3-RelativeResponse.xlsx"#filter functions
#Input and output files
read.inname = "./OutData/AOT/input_file.txt"#where to save/call input file
read.outname = "./OutData/AOT/output_" #where to save/call output file
#Use this file to define your AOT
out.sum.table = "./OutData/AOT/summaryoutput.txt"#This file will show you how the reflectance of the optically deepwater changes with AOT

#Define your ^SVinput file
tab = rep("",17)
tab[1] = "0" #User Defined
#These are image specific parameters that must be pulled from the metadata
tab[2] = "33.9 152.5 28.1 133.6 8 17"# (solar zenithal, solar azimuthal, sensor zenithal, senzor azimuthal, month, day)
#
tab[3] = "2" #(Midlatitude summer)
tab[4] = "2" #(Maritime Model)
tab[5] = "0"
tab[6] = "0.10" #value (aot)
tab[7] = "0" #(target level)
tab[8] = "-1000" #(sensor level, satellite)
tab[9] = "1" #(User's defined filtered function)
#Read in new filter function
excel.in = read_xlsx(dir.relative.response,sheet=band.num)#read in excel
tab[10] = paste(as.numeric(excel.in[1,1])/1000,as.numeric(excel.in[dim(as.matrix(excel.in[,1]))[1],1])/1000,collapse=" ")#(start and stop wavelengths)
tab[11] = paste(round(as.matrix(excel.in[,2]),3),collapse=" ")#Filter Function
tab[12] = "0" #Homogeneous surface
tab[13] = "1" #(directional effects)
tab[14] = "6" #Ocean
tab[15] = "3.89 220 34 0.5"#(Wind speed (m/s) Wind Azim. (in degrees) salinity (degrees) pigment concentration (mg/m3))
tab[16] = "1" #BRDF
tab[17] = "-0.0129" #min reflectance of deepwater pixel 
##
rm(dir.relative.response,excel.in)

##Define your output summary table
a=matrix(NA, ncol=7, nrow=1)
colnames(a) = c("AOT", "Min", "1st", "Median","Mean","3rd","Max")
write.table(a, file=out.sum.table,row.names = F)
rm(a)
#

##Loop through
aot = as.character(seq(0.001,0.22,0.001))#What AOT to use
for (i in aot){
  print(i)
  tab[6] = i #change aot
  #Run 6SV
  write.table(tab,read.inname, row.names = F, col.names = F, quote=F)#write new input file
  cmd = paste(paste("sixsV2.1 <",read.inname,sep=""), paste(">",read.outname,i,".txt",sep=""))#System command
  shell(cmd)#Run 6SV
  #Extract coefficients
  sixSV = read.delim(paste(read.outname,i,".txt",sep=""),header=F,colClasses = "character",fill=T,sep=c(""),skipNul = T,na.strings = "")
  a.ref = which(sixSV == "xap", arr.ind=TRUE)#find position of ref coefficients 
  print(sixSV[c(as.numeric(a.ref[1]):as.numeric(a.ref[1]+1)),])
  a.ref=sixSV[a.ref[1]+1,]#extract row of data
  a.ref = a.ref[, !apply(is.na(a.ref), 2, all)]#remove NAS
  if (a.ref[1]!=":"){
  xap = as.numeric(a.ref[1])
  xb = as.numeric(a.ref[2])
  xc = as.numeric(a.ref[3])}
  if (a.ref[1]==":"){
    xap = as.numeric(a.ref[2])
    xb = as.numeric(a.ref[3])
    xc = as.numeric(a.ref[4])}
  #Calculate BOA reflectance
  y.ref=(xap*nir.ref)-xb #6SV formula
  acr.ref=y.ref/(1.+(xc*y.ref))#6SV formula
  #Summary stats
  sum.acr.ref = t(do.call(cbind, lapply(acr.ref, summary)))
  row.names(sum.acr.ref) = i
  write.table(sum.acr.ref, file=out.sum.table,append=T,col.names = F) 
  rm(acr.ref,cmd,sum.acr.ref,xap,xb,xc,y.ref,sixSV,a.ref)}






