load("~/LAVA/data/x.RData")

param <- commandArgs(trailingOnly=T)
i= eval(paste(text=param[1]))
dir.create(i)
i=as.numeric(i)
print(i)

library(data.table)
library(R.utils)

chr=loci[i,2]
start=loci[i,3]
stop=loci[i,4]

plinkpref="~/Plink --bfile ~/ldsc/partition/1000G.EUR.QC."
plinksuf=" --make-bed --noweb --out"
plinkcode=paste(c(plinkpref,chr," --chr ",chr," --from-bp ",start," --to-bp ",stop,plinksuf," ",i,"/ref"),collapse="")
print(plinkcode)
system(plinkcode)
ref=fread(paste(i,"ref.bim",sep="/"),data.table=F)
ref=nrow(ref)
int=loci_heritability[which(loci_heritability$id==i),]
pheno=int[which(int$p<3.17e-5),1]

filepref="~/ldsc/all_gwas/"
for(dis in pheno){
  file=paste(c(filepref,dis,".gz"),collapse="")
  d=fread(file,data.table=F,fill=TRUE)
  d=d[which(d$CHR==chr),]
  d=d[which(d$BP>start),]
  d=d[which(d$BP<stop),]
  d$Z=d$b/d$se
  d=d[,c("SNP","A1","A2","CHR","BP","Z","N")]
  d=d[which(!is.na(d$Z) & is.finite(d$Z)),]
  if(nrow(d)>0.3*ref){fwrite(d,paste(c(i,"/",dis,".gz"),collapse=""),sep=" ")}
}
system(paste(c("./fiziscript",chr,start,stop,i),collapse=" "))

  