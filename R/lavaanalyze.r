library(data.table)
disname=fread("disname.csv",data.table=F)
lava=fread("res_lava.csv",data.table=F)
lava$name1=disname[lava[,1],2]
lava$name2=disname[lava[,2],2]
lava$fdr=p.adjust(lava$p/2,method="fdr")
loci_antag=lava[which(lava$rg<0&lava$p<5e-7),]
loci_antag=lava[which(lava$rg<0&lava$fdr<0.05),]
head(loci_antag)

length(unique(loci_antag$id))y67
length(union(loci_antag[,1],loci_antag[,2]))
length(which(loci_antag[,1] %in% red & loci_antag[,2]
antag_loci_bed=list()

for(i in disname[,1]){
	index=which(loci_antag[,1]==i | loci_antag[,2]==i)
	if(length(index)==0){next}
	id=loci_antag[index,"id"]
	bed=loci[id,2:4]
	antag_loci_bed[[i]]=bed
}
names(antag_loci_bed)

all_loci_bed=list()

for(i in disname[,1]){
	index=which(loci_heritability[,1]==i&loci_heritability$p<3.17e-5)
	if(length(index)==0){next}
	id=loci_heritability[index,"id"]
	bed=loci[id,2:4]
	all_loci_bed[[i]]=bed
}


pref="F:\\selection\\newbaseline\\baselineLD."
for(i in 1:22){
	file=paste(c(pref,i,".annot.gz"),collapse="")
	d=fread(file,data.table=F)
	d=d[,c("BP","SNP","Recomb_Rate_10kb")]
	rec[[i]]$REC=d[match(rec[[i]]$SNP,d$SNP),3]
	#rec[[i]][which(rec[[i]]$CM<0),"CM"]=0
	#rec[[i]]$CM=floor(rec[[i]]$CM*100)
}

twoway=data.frame(rbindlist(rec))
thred=list()
for(metric in c("B2","ASMC","GERP","CMS","SDS","GROSS")){
	thred[[metric]]=quantile(twoway[which(is.finite(twoway[,metric])),metric],0.99)
}	


loci$rec=NA
loci$lld=NA
loci$div=NA
loci$age=NA
loci$GREP=NA
loci$ASMC=NA
loci$CMS=NA
loci$SDS=NA
loci$GROSS=NA



for(i in 1:nrow(loci)){
	chr=loci[i,2]
	start=loci[i,3]
	end=loci[i,4]
	int=rec[[chr]]
	int=int[which(int$BP>start),]
	int=int[which(int$BP<end),]
	int1=int[which(int[,3]!=0),3]
	loci[i,"rec"]=median(int1)
	int1=int[which(int[,"AGE"]!=0),"AGE"]
	loci[i,"age"]=median(int1)
	int1=int[which(int[,"DIV"]!=0),"DIV"]
	loci[i,"div"]=median(int1)
	int1=int[which(int[,"LLD"]!=0),"LLD"]
	loci[i,"lld"]=median(int1)
	int1=int[which(int[,"GERP"]!=0),"GERP"]
	loci[i,"GERP"]=max(int1)
	int1=int[which(int[,"ASMC"]!=0),"ASMC"]
	loci[i,"ASMC"]=max(int1)
	int1=int[which(int[,"cms"]!=0),"cms"]
	loci[i,"CMS"]=max(int1)
	int1=int[which(int[,"sds"]!=0),"sds"]
	loci[i,"SDS"]=max(int1)	
	int1=int[which(int[,"gross"]!=0),"gross"]
	loci[i,"GROSS"]=max(int1)	
}

res_antag_loci=c()
for(dis in disname[,1]){
	index=which(loci_heritability[,1]==dis&loci_heritability$p<3.17e-5)
	al=loci_heritability[index,"id"]
	index=which(loci_antag[,1]==dis | loci_antag[,2]==dis)
	tl=loci_antag[index,"id"]	
	res=sapply(c("GERP","ASMC","CMS","SDS","GROSS"),
		function(x){get_tva_t(al=al,tl=tl,metric=x)})
	res_antag_loci=rbind(res_antag_loci,c(dis,res))
}

for(i in 1:22){
	x=rec[[i]]
	chr=x[1,1]
	int=loci[which(loci[,2]==chr),]	
	breaks=c(1,int[,4])
	labels=int[,1]
	res <- x %>% mutate(loci=cut(BP, breaks=breaks, labels=labels))
	rec[[i]]=res
}

binary_loci_cm=c()
for(dis in disname[5:111,1]){
	index=which(loci_heritability[,1]==dis&loci_heritability$p<3.17e-5)
	al=loci_heritability[index,"id"]
	index=which(loci_antag[,1]==dis | loci_antag[,2]==dis)
	tl=loci_antag[index,"id"]	
	if(length(tl)<5){next}
	int=c()
	for(i in 1:22){
		int=rbind(int,rec[[i]][which(rec[[i]]$loci %in% al),])
	}
	int$CM=paste(int$CHR,int$CM,sep="_")
	int$group=ifelse(int$loci %in% tl,"yes","no")
	int$loci=as.numeric(int$loci)
	int1=int %>% group_by(CM) %>% summarise(
		rec=median(REC,na.rm=T),
		age=median(AGE,na.rm=T),
		lld=median(LLD,na.rm=T),
		div=median(DIV,na.rm=T),
		ASMC=max(ASMC,na.rm=T),
		GERP=max(GERP,na.rm=T),
		CMS=max(cms,na.rm=T),
		SDS=max(sds,na.rm=T),
		GROSS=max(gross,na.rm=T),
		B2=max(B2,na.rm=T),
		n=n(),
		loci=median(loci)) %>% arrange(CM)
	int1=data.frame(int1)
	int1$group=ifelse(int1$loci %in% tl,"yes","no")
	res=c()
	for(metric in c("ASMC","GERP","SDS","CMS","GROSS","B2")){
		dmod=int1[which(is.finite(int1[,metric]) & !is.na(int1[,metric])),
				c("rec","age","lld","div","n",metric,"group")]
		dmod[,6]=ifelse(dmod[,6]>thred[[metric]],1,0)
		if(nrow(dmod)==0){next}
		colnames(dmod)[6]="metric"
		model=glm(metric~.,family=binomial(),data=dmod)
		z=summary(model)$coefficients["groupyes",3]
		res=rbind(res,c(metric,z))
	}
	res=data.frame(res,dis=dis)
	binary_loci_cm=rbind(binary_loci_cm,res)
}

#############################################
#################PRS#########################
#############################################
library(ggplot2)
prs=fread("prs.txt",data.table=F,header=F)
colnames(prs)=c("dis","study","type","source","beta","p")

getD=function(study,source){
	d=prs[which(prs$study==study & prs$source==source),]
	d$beta=d$beta/abs(d$beta)
	d[which(d$p==0),"p"]=0.001
	d$p=abs(qnorm(d$p))
	d$p=d$p*d$beta
	colnames(d)=c("dis","study","type","source","sign","z")
	return(d)
}

draw=getD("loci","relate")
ggplot(draw,aes(z,type,group=type,fill=type))+geom_boxplot()		
prs$name=disname[match(prs[,1],disname[,1]),2]
write.csv(prs,"prs.csv")




