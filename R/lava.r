load("~/LAVA/data/x.RData")
param <- commandArgs(trailingOnly=T)
i= eval(paste(text=param[1]))
setwd(i)

library(data.table)
library(R.utils)
library(LAVA)
system("rm *sumstat")
system("rm *log")
tmp=list.files(pattern="*gz")

tmp=sub(".gz","",tmp)
input = process.input(input.info.file="../info.txt", 
  ref.prefix="ref",
  sample.overlap.file="../sample.overlap.txt",
  phenos=tmp) 
locus = process.locus(loci[i,], input)
print(i)
res=run.bivar(locus)
write.table(res,"lava.res",quote=F,row.names=F)