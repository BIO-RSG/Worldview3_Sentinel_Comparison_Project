##
rm(list=ls())
library(readxl)
#
setwd("D:\\SixSV")#For R
shell("cd D:\\SixSV")#For command line script
#
dir.relative.response = "./OutData/WV3-RelativeResponse.xlsx"#filter functions
read.inname = "./OutData/20190811/AC/input_file_"#where to save/call input file
read.outname = "./OutData/20190811/AC/output_" #where to save/call output file
#

for (i in 1:8){
  #Define your input file
  tab = rep("",17)
  tab[1] = "0" #User Defined
  tab[2] = "33.9 152.5 28.1 133.6 8 17"# (solar zenithal, solar azimuthal, sensor zenithal, senzor azimuthal, month, day)
  tab[3] = "2" #(Midlatitude summer)
  tab[4] = "2" #(Maritime Model)
  tab[5] = "0"
  tab[6] = "0.08" #value (aot)
  tab[7] = "0" #(target level)
  tab[8] = "-1000" #(sensor level, satellite)
  tab[9] = "1" #(User's defined filtered function)
  #Read in new filter function
  excel.in = read_xlsx(dir.relative.response,sheet=i)#read in excel
  tab[10] = paste(as.numeric(excel.in[1,1])/1000,as.numeric(excel.in[dim(as.matrix(excel.in[,1]))[1],1])/1000,collapse=" ")#(start and stop wavelengths)
  tab[11] = paste(round(as.matrix(excel.in[,2]),3),collapse=" ")#Filter Function
  tab[12] = "0" #Homogeneous surface
  tab[13] = "1" #(directional effects)
  tab[14] = "6" #Ocean
  tab[15] = "3.89 220 34 0.5"#(Wind speed (m/s) Wind Azim. (in degrees) salinity (degrees) pigment concentration (mg/m3))
  tab[16] = "1" #BRDF
  tab[17] = "-0.1" #min reflectance of deepwater pixel 
  ##
  #Run 6SV
  write.table(tab,paste(read.inname,i,".txt",sep=""), row.names = F, col.names = F, quote=F)#write new input file
  cmd = paste(paste("sixsV2.1 <",read.inname,i,".txt",sep=""), paste(">",read.outname,i,".txt",sep=""))#System command
  shell(cmd)#Run 6SV
  }






