ibrary(data.table)
library(R.utils)
library(readr)
memory.limit(999999)
library(tidyr)

N=data.frame(N,thred=0.05/as.numeric(N[,2]))
head(N)
rownames(N)=N[,1]
allpred=data.frame(allpred,thred=N[allpred$name,3])
head(allpred)
allpred=allpred[which(allpred$pvalue<allpred$thred),]
int=t(data.frame(strsplit(allpred$name,"+",fixed=T)))
allpred=data.frame(allpred,trait=int[,1],tissue=int[,2])
allpred=allpred[,-7]
allpred[,1]=t(data.frame(strsplit(allpred[,1],".",fixed=T)))[,1]

############################################################
N=data.frame(N,thred=0.05/as.numeric(N[,2]))
head(N)
rownames(N)=N[,1]

allsmr=data.frame(allsmr,thred=N[allsmr$name,3])
head(allsmr)
allsmr=allsmr[which(allsmr[,4]<allsmr$thred),]
int=t(data.frame(strsplit(allsmr$name,"+",fixed=T)))
allsmr=data.frame(allsmr,trait=int[,2],tissue=int[,1])
allsmr=allsmr[,-6]
allsmr[,1]=t(data.frame(strsplit(allsmr[,1],".",fixed=T)))[,1]
l=which(is.na(allsmr[,2]))
allsmr[l,2]=allsmr[l,1]
allsmr=allsmr[which(allsmr[,3]!="inf"&allsmr[,3]!="-inf"),]
all_smr=allsmr[which(all_smr[,3]!="inf"&all_smr[,3]!="-inf"),]

focusres=focusres[!is.na(focusres$tissue),]
###############################################################
antagres=data.frame(gene=c(),T1=c(),b1=c(),t1=c(),e1=c(),T2=c(),b2=c(),t2=c(),e2=c())
all_smr=allsmr[which(allsmr[,5]>0.05),]

for(T in disname[,1]){
	int=get_assoc(T)
	int=data.frame(format_int(int))
	int=int[which(int$b1!=0),]
	intT=get_antag(int)
	antagres=rbind(antagres,intT)
}

save(antagres,file="antagres.RData")

intres=antagres[!duplicated(paste(antagres[,1],antagres[,2],sep="")),]
table(intres$T1)

all_assoc=c()
for(T in disname[,1]){
	int=get_assoc(T)
	int=data.frame(format_int(int))
	int=int[which(int$b1!=0),]
	all_assoc=rbind(all_assoc,int)
}
x=unique(all_assoc[grep("ENSG",all_assoc[,1]),1])
gene.d <-bitr(x, fromType = "ENSEMBL", 
        toType = "SYMBOL",
        OrgDb = org.Hs.eg.db)
xx=setdiff(x,gene.d[,1])
gene.d=rbind(gene.d,data.frame(ENSEMBL=xx,SYMBOL=xx))
gene.d=gene.d[!duplicated(gene.d[,1]),]
rownames(gene.d)=gene.d[,1]
all_assoc[grep("ENSG",all_assoc[,1]),1]=gene.d[all_assoc[grep("ENSG",all_assoc[,1]),1],2]
antagres[grep("ENSG",antagres[,1]),1]=gene.d[antagres[grep("ENSG",antagres[,1]),1],2]


table(as.vector(table(all_assoc[,1])))

disname=fread("disname.csv",data.table=F)
rownames(disname)=disname[,1]

l=table(all_assoc[,1])
table(as.vector(l))
L=l[which(!names(l) %in% mhcgene)]
l[which(l>44)]
L[which(L>30)]
disname[all_assoc[which(all_assoc[,1]=="C4A"),2],2]

disname["ID636",2]="SLE"

gene_stat=c()
for(gene in unique(all_assoc[,1])){
n=length(which(all_assoc[,1]==gene))
neg=length(which(all_assoc[,1]==gene&all_assoc[,3]<0))
mhc=ifelse(gene %in% mhcgene,1,0)
gene_stat=rbind(gene_stat,c(gene=gene,n=n,neg=neg,mhc=mhc))
}
gene_stat=data.frame(gene_stat)
gene_stat[,2:4]=apply(gene_stat[,2:4],2,as.numeric)
gene_stat=gene_stat[order(gene_stat[,2],decreasing=TRUE),]
head(gene_stat[which(gene_stat[,4]==0),],n=50)

gene_stat$loeuf=loeuf[match(gene_stat[,1],loeuf[,1]),2]
gene_stat$hi=HI[gene_stat[,1],5]
cons=loeuf[which(loeuf[,2]<0.35),1]
spearman_test(loeuf~n,data=
gene_stat[which(gene_stat$mhc==0&!is.na(gene_stat$loeuf)),])

dis_stat=c()
for(dis in unique(antagres$T1)){
gene=paste(unique(antagres[which(antagres$T1==dis|antagres$T2==dis),1]),collapse=" ")
n=length(unique(antagres[which(antagres$T1==dis),1]))
dis_stat=rbind(dis_stat,c(dis=dis,n=n,gene=gene))
}
dis_stat=data.frame(dis_stat)
dis_stat=dis_stat[order(as.numeric(dis_stat$n),decreasing=TRUE),]
dis_stat$n=as.numeric(dis_stat$n)
dis_stat$name=disname[match(dis_stat[,1],disname[,1]),2]
dis_stat=dis_stat[!is.na(dis_stat[,2]),]
cons=loeuf[which(loeuf[,2]<0.35),1]
hi=rownames(HI)[18000:18263]

a=gene_stat[which(gene_stat[,1] %in% cons),2:3]
t.test(a[,2]/a[,1])

unique(antagres[,2])
for(dis in unique(antagres[,2])){
gene=paste(unique(antagres[which(antagres[,2]==dis),1]),collapse=" ")
write.table(data.frame("antag",gene),paste("./disantag/",dis,sep=""),quote=F,row.names=F,col.names=F)
}

antaggene=unique(antagres[,1])
save(antaggene,file="antaggene.RData")

disease=disname[,1]
dis_network=c()
for(i in 1:(length(disease)-1)){
for(j in (i+1):length(disease)){
T1=disease[i]
T2=disease[j]
N=intersect(unique(all_assoc[which(all_assoc[,2]==T1),1]),
unique(all_assoc[which(all_assoc[,2]==T2),1]))
Nnonmhc=setdiff(N,mhcgene)
N=length(N)
Nnonmhc=length(Nnonmhc)

n1=antagres[which(antagres$T1==T1&antagres$T2==T2),1]
n2=antagres[which(antagres$T1==T2&antagres$T2==T1),1]
n=unique(c(n1,n2))
nnonmhc=setdiff(n,mhcgene)
n=length(n)
nnonmhc=length(nnonmhc)
dis_network=rbind(dis_network,c(T1=T1,T2=T2,N=N,n=n,Nnonmhc=Nnonmhc,nnonmhc=nnonmhc))
}
}

dis_network=data.frame(dis_network)
dis_network[,3:6]=apply(dis_network[,3:6],2,as.numeric)
head(dis_network[order(dis_network[,6],decreasing=TRUE),],n=20)
disname["ID3",]

a=as.numeric(all_assoc[which(all_assoc[,1] %in% hi),3])
b=a[which(a<0)]
binom.test(length(b),length(a))

for(i in 1:nrow(gs)){
	gene=unique(all_assoc[which(all_assoc[,2]==gs[i,1]),1])
	n=length(gene)
	gene1=paste(gene,collapse=" ")
	gs[i,3]=gene1

	gene2=unique(antagres[which(antagres[,2]==gs[i,1]),1])
	gene2=intersect(gene,gene2)
	N=length(gene2)
	gene2=paste(gene2,collapse=" ")
	gs[i,4]=gene2
	
	gs[i,2]=paste(N,n,sep="/")
}

#all assoc gene, all pleiotropic gene, all antaggene, all antag- and reversly adapted gene

genelist=c()
for(id in rownames(disname)){
	gene=unique(all_assoc[which(all_assoc[,2]==id),1])
	n1=length(gene)
	gene1=paste(gene,collapse=" ")

	gene2=intersect(gene,plegene)
	n2=length(gene2)
	gene2=paste(gene2,collapse=" ")

	g3=intersect(gene,antaggene)
	n3=length(g3)
	gene3=paste(g3,collapse=" ")	

	g3=intersect(g3,rownames(sum_traj_500))
	int=all_assoc[which(all_assoc[,2]==id),]
	agl=as.numeric(int[match(g3,int[,1]),3])
	agl1=sign(sum_traj_500[g3,"s500"])
	gene4=rownames(sum_traj_500)[which(sum_traj_500[g3,"p500"]<0.05&tgl*tgl1<0)]
	n4=length(gene4)
	gene4=paste(gene4,collapse=" ")

	genelist=rbind(genelist,c(id,disname[id,2],gene1,n1,gene2,n2,gene3,n3,gene4,n4))
}
write.csv(genelist,"geneist.csv")
