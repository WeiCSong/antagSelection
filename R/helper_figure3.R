get_top_signal=function(gene,metric){
	d=twoway
	if (metric=="gross") {d=twgross}
	snp=g2q[which(g2q[,2]==gene),1]
	int=d[which(d$SNP %in% snp),metric]
	int=max(int)
	return(int)
}

get_tva_t=function(ag,tg,metric){
	t=0
	if(length(ag)>length(tg) & length(tg)>10){
		d=top_signal_gene[which(rownames(top_signal_gene) %in% ag),]
		d=d[which(d[,metric]!=0 & is.finite(d[,metric])),]
		d=data.frame(d[,c("AGE","LLD","DIV","nsnp")],
				 metric=d[,metric],
				 group=ifelse(rownames(d) %in% tg,"yes","no"))
		model=lm(metric~.,data=d)
		t=summary(model)$coefficients["groupyes",3]
	}
	return(c(metric,t))
}

#	
metric2gene=function(metric){
	d=twoway
	if(metric=="gross"){d=twgross}
	d=d[which(d[,metric]!=0),c(metric,"SNP")]
	threshold=quantile(d[,1],0.99)
	snp=d[which(d[,1]>threshold),2]
	snp1=intersect(snp,g2q[,1])
	snp2=setdiff(snp,snp1)
	qgene=unique(g2q[which(g2q[,1] %in% snp1),2])
	snp.info=SNP[which(SNP$Name %in% snp2),]
	res=snps.to.genes(snp.info, gene, 100000)	
	pgene=names(which(sapply(res,length)!=0))
	gene=unique(c(pgene,qgene))
	return(gene)
}
#

snp2avgMetric=function(snp,metric){
	d=twoway
	if(metric=="gross"){d=twgross}
	index=which(d$SNP==snp)
	bp=d[index,"BP"]
	value=c()
	index1=index
	while(d[index1,"BP"]+100000>bp){
		value=c(value,d[index1,metric])
		index1=index1-1
	}
	index1=index
	while(d[index1,"BP"]-100000<bp){
		value=c(value,d[index1,metric])
		index1=index1+1
	}
	value=median(value[which(value!=0)])
}

#
snp2avgMetric=function(snp,metric){
	d=twoway
	if(metric=="gross"){d=twgross}
	index=which(d$SNP==snp)
	d=d[which(d$CHR==d[index,"CHR"]),]
	d=d[which(d$BP>d[index,"BP"]-100000),]
	d=d[which(d$BP<d[index,"BP"]+100000),]	
	d=d[which(d[,metric]!=0),]
	value=median(d[,metric])
}
#



top_metric=c()
for(i in 1:22){	
	int=rec[[i]]

	int1=int %>% group_by(CM) %>% summarise(
		B2=max(B2,na.rm=T),
		ASMC=max(ASMC,na.rm=T),
		GERP=max(GERP,na.rm=T),
		CMS=max(cms,na.rm=T),
		SDS=max(sds,na.rm=T),
		GROSS=max(gross,na.rm=T),
		start=min(BP),
		end=max(BP)) %>% arrange(CM)
	int1=data.frame(int1)
	res=c()
	for(metric in c("B2","ASMC","GERP","CMS","SDS","GROSS")){
		int=int1[which(!is.na(int1[,metric]) & is.finite (int1[,metric])),]
		int=int[which(int[,metric]>thred[[metric]]),]
		rownames(int)=as.character(int$CM)
		bk=int
		
		for(rn in rownames(int)[-1]){
			index=match(rn,rownames(int))
			index1=match(rn,rownames(bk))
			bk[index1-1,"end"]=bk[index1-1,"end"]+1000
			bk[index1,"start"]=bk[index1,"start"]-1000
			if(int[index,"CM"]<bk[index1-1,"CM"]+20){
				bk[index1-1,"end"]=bk[index1,"end"]	
				bk[index1-1,"CM"]=bk[index1,"CM"]	
				bk=bk[-index1,]
			}
		}
		bk=data.frame(metric=metric,start=bk$start,end=bk$end)
		res=rbind(res,bk)
	}
	res=data.frame(chr=i,res)
	top_metric=rbind(top_metric,res)
}			

metric_regression=function(int1=int1,method=c("continuous","binary")){
	res=c()
	for(metric in c("ASMC","GERP","SDS","CMS","GROSS","B2")){
		dmod=int1[which(is.finite(int1[,metric]) & !is.na(int1[,metric])),
				c("rec","age","lld","div","n",metric,"group")]
		if(nrow(dmod)==0){next}
		if(method=="binary"){
			dmod[,6]=ifelse(dmod[,6]>thred[[metric]],1,0)
			colnames(dmod)[6]="metric"
			model=glm(metric~.,family=binomial(),data=dmod)
			t=summary(model)$coefficients["groupyes",3]
		}else{
			colnames(dmod)[6]="metric"
			model=lm(metric~.,data=dmod)
			t=summary(model)$coefficients["groupyes",3]	
		}		
		res=rbind(res,c(metric,t))
	}
	return(data.frame(res))
}




