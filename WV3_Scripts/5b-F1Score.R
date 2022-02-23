
#Calculate F1 score
fun.f1 = function(tp, fp, fn){
  fun.pr = tp/(tp+fp)
  fun.rc = tp/(tp+fn)
  (2*fun.pr*fun.rc)/(fun.pr+fun.rc)
  
}#F1
fun.tp = function(in.csv){
  dat.in = paste0(in.csv,list.files(in.csv, pattern=".csv$"))
  b = rep(NA, length(dat.in))
  for (i in dat.in){
    a = read.csv(i)
    b[which(dat.in==i)]=a[1,1]
  }
  b
}#True positive
fun.fp = function(in.csv){
  dat.in = paste0(in.csv,list.files(in.csv, pattern=".csv$"))
  b = rep(NA, length(dat.in))
  for (i in dat.in){
    a = read.csv(i)
    b[which(dat.in==i)]=a[2,1]
  }
  b
}#False negative
fun.fn = function(in.csv){
  dat.in = paste0(in.csv,list.files(in.csv, pattern=".csv$"))
  b = rep(NA, length(dat.in))
  for (i in dat.in){
    a = read.csv(i)
    b[which(dat.in==i)]=a[1,2]
  }
  b
}##False negative
##


#Aug 17
in.dir = "C:\\Where\\is\\my\\data\\"
dir.use = list.dirs(in.dir, recursive=F)
dir.use = dir.use[c(7,6,5,14,13,12,11,8,15,1)]
input.lab = c("6B","5B","4B","PC1-5","PC1-4","PC1-3","PC1-2","DII", "6B+PC1-4+DII", "B-G-R")
out.f1 = matrix(NA, ncol=length(input.lab), nrow=50)
colnames(out.f1) = input.lab
for (j in dir.use){
  print(j)
  in.csv = paste0(j,"//cvrunmf//")
  tp = fun.tp(in.csv)
  fp =fun.fp(in.csv)
  fn = fun.fn(in.csv)
  f1 = fun.f1(tp,fp,fn)
  out.f1[,which(dir.use==j)] = f1
}
aug.17.f1 = apply(out.f1,2, mean)

#Aug 11
in.dir = "C:\\Where\\is\\my\\data\\"
dir.use = list.dirs(in.dir, recursive=F)
dir.use = dir.use[c(7,6,5,15,12,11,9,8,17,18,10,1)]
input.lab = c("6B","5B","4B","PC1-5","PC1-4","PC1-3","PC1-2","DII","6B+PC1-4+DII","6B+PC1-4+DII+SDB","BG+PC12+SDB","B-G-R")
out.f1 = matrix(NA, ncol=length(input.lab), nrow=50)
colnames(out.f1) = input.lab
for (j in dir.use){
  print(j)
  in.csv = paste0(j,"//cvrunmf//")
  tp = fun.tp(in.csv)
  fp =fun.fp(in.csv)
  fn = fun.fn(in.csv)
  f1 = fun.f1(tp,fp,fn)
  out.f1[,which(dir.use==j)] = f1
}
aug.11.f1 = apply(out.f1,2, mean)

#S-2
in.csv="C:\\Where\\is\\my\\data\\RFResultv2//FinalModels//S2-v2//cvrun//"
tp = fun.tp(in.csv)
fp =fun.fp(in.csv)
fn = fun.fn(in.csv)
f1 = fun.f1(tp,fp,fn)
s2.f1 = mean(f1)
