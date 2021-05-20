#This is an R implement of focus import function. It aviod the [ERROR "symbol"] and the query of gencode database, thus is computationally efficient (~0.6h for importing all gtexv8 elastic net model).
#Since the tx_start and tx_end of genes does not actually matter, i replaced them by the midpoint of correspponding eqtls, which ensure that all cis-eqtl will be included. 

library(tidyr)
setwd("/lustre/home/acct-bmelgn/bmelgn-3/focusdata/")
library(RSQLite)
load("0519.RData") #contains: df kgsnp, which match rsid to hg19 chr and pos; empty df refpanel, model, model attribute, molecularfeature and weight, with the same colnames as focus.db provided by focus author.
tmp=list.files(pattern="en*")
###read in .db
for (filename in tmp){
tissue=sub("en_","",filename)
tissue=sub(".db","",tissue)

sqlite.driver <- dbDriver("SQLite")
db1 <- dbConnect(sqlite.driver,
                dbname = filename)
weights=dbReadTable(db1,"weights")
extra=dbReadTable(db1,"extra")

weights=weights[which(weights[,2] %in% kgsnp[,2]),] #I restrict my analysis to hg19 database, due to the error of 1000G liftover (see their ftp site)
l=match(weights[,2],kgsnp[,2])

weights$chr=paste("chr",kgsnp[l,1],sep="")
weights$pos=kgsnp[l,4]
weights[,1]=t(data.frame(strsplit(weights[,1],".",fixed=T)))[,1]

############################update refpanel######################################
refpanel=rbind(refpanel,data.frame(id=max(refpanel[,1])+1,ref_name="gtex",tissue=tissue,assay="rnaseq"))
rownames(refpanel)=1:nrow(refpanel)

###########################update molecularfeature#################################
int2=data.frame(id="id",ens_gene_id=t(data.frame(strsplit(extra[,1],".",fixed=T)))[,1],ens_tx_id=NA,mol_name=extra[,2],type=extra[,3])
gene=int2[,2]
rownames(int2)=c()
int2=int2[which(! int2[,2] %in% molecularfeature[,2]),]
f=function(gene){
x=weights[which(weights[,1]==gene),"pos"]
return(mean(c(min(x),max(x))))
}

int2=data.frame(int2,chrom=weights[match(int2[,2],weights[,1]),"chr"],tx_start=unlist(sapply(int2[,2],f)))
int2$tx_start=ceiling(int2$tx_start)
int2=data.frame(int2,tx_end=int2$tx_start+1)

molecularfeature=rbind(molecularfeature,int2)
molecularfeature$id=1:nrow(molecularfeature)
rownames(molecularfeature)=1:nrow(molecularfeature)

########################################update model#####################################
int3=data.frame(id="id",inference="en",ref_id=nrow(refpanel),mol_id=match(gene,molecularfeature[,2]))
model=rbind(model,int3)
model$id=1:nrow(model)
rownames(model)=1:nrow(model)
intmodel=model[which(model[,3]==nrow(refpanel)),]

#############################update modelattribute######################################
extra[,1]=t(data.frame(strsplit(extra[,1],".",fixed=T)))[,1]
extra[,1]=molecularfeature[match(extra[,1],molecularfeature[,2]),1]
extra[,1]=intmodel[match(extra[,1],intmodel[,4]),1]
intextra=extra[,c(1,9,12)]
colnames(intextra)=c("model_id","cv.R2","cv.R2.pval") #if using mashr model, please modify the column name.
int4=gather(intextra,attr_name,value,cv.R2:cv.R2.pval)
int4=int4[order(int4[,1]),]
int4=data.frame(id="id",int4[,c(2,3,1)])
modelattribute=rbind(modelattribute,int4)
modelattribute$id=1:nrow(modelattribute)
rownames(modelattribute)=1:nrow(modelattribute)

################################update weight##############################
int5=weights[,c(2,7,8,5,4,6)]
int5=data.frame(int5,se=NA)
int5$model=molecularfeature[match(weights[,1],molecularfeature[,2]),1]
int5$model=intmodel[match(int5$model,intmodel[,4]),1]
int5=data.frame(id="id",int5)
colnames(int5)=c("id","snp","chrom","pos","effect_allele","alt_allele","effect","std_error","model_id")
weight=rbind(weight,int5)
weight$id=1:nrow(weight)
rownames(weight)=1:nrow(weight)
}

rm(kgsnp)
filename <- "gtex_v8.db"
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver,
                dbname = filename)
dbWriteTable(db,"refpanel",refpanel,overwrite=T)
dbWriteTable(db,"molecularfeature",molecularfeature,overwrite=T)
dbWriteTable(db,"model",model,overwrite=T)
dbWriteTable(db,"modelattribute",modelattribute,overwrite=T)
dbWriteTable(db,"weight",weight,overwrite=T)
