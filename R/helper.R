sign=function(x){
	x/abs(x)
}

get_assoc=function(T){
	int1=allsmr[which(allsmr$trait==T),c(1,2,7,3,8,4)]
	colnames(int1)=c("id","gene","T1","b1","t1","p")
	if(nrow(int1)!=0){
		int1=data.frame(int1,e1="SMR")
	}

	int2=allpred[which(allpred$trait==T),c(1,2,8,3,9,4)]
	colnames(int2)=c("id","gene","T1","b1","t1","p")
	if(nrow(int2)!=0){
		int2=data.frame(int2,e1="pred")	
	}	

	int3=focusres[which(focusres$trait==T),c(1,2,6,4,3,5)]
	int3[,6]=1-int3[,6]
	colnames(int3)=c("id","gene","T1","b1","t1","p")
	if(nrow(int3)!=0){	
		int3=data.frame(int3,e1="focus")
	}

	int=rbind(int1,int2)
	int=rbind(int,int3)
	if(nrow(int)>0){
		int[,c(4,6)]=apply(int[,c(4,6)],2,as.numeric)
	}
	return(int)
}



format_int=function(int){
	out=c()
	for (gene in unique(int[,1])){
		INT=int[which(int[,1]==gene),]
		if(max(na.omit(sign(INT$b1)))!=min(na.omit(sign(INT$b1)))){
			next
		}
		e1=paste(unique(INT[,"e1"]),collapse="|")
		x=which(INT$p==min(INT$p))[1]
		b1=INT[x,"b1"]
		t1=paste(unique(INT[,"t1"]),collapse="|")
		T1=INT[x,"T1"]
		out=rbind(out,c(gene=gene,T1=T1,b1=b1,t1=t1,e1=e1))
	}
	return(out)
}

get_antag=function(int){
	out=c()
	for(x in 1:nrow(int)){
		gene=int[x,1]
		T1=int[x,2]
		s=sign(as.numeric(int[x,3]))
		avadis=c(unique(all_smr[which(all_smr[,1]==gene),"trait"]),unique(allpred[which(allpred[,1]==gene),"trait"]))
		avadis=intersect(setdiff(unique(c(avadis,unique(focusres[which(focusres[,1]==gene),"trait"]))),T1),dis)
		if(length(avadis)==0){
			next
		}
		intout=c()
		for(T2 in avadis){
			index=which(all_smr$trait==T2)
			int1=all_smr[index,c(1,7,3,8,4)]
			colnames(int1)=c("gene","T1","b1","t1","p")
			if(nrow(int1)!=0){
				int1=data.frame(int1,e1="SMR")
			}

			int2=allpred[which(allpred$trait==T2),c(1,8,3,9,4)]
			colnames(int2)=c("gene","T1","b1","t1","p")
			if(nrow(int2)!=0){
				int2=data.frame(int2,e1="pred")	
			}

			int3=focusres[which(focusres$trait==T2),c(1,6,4,3,5)]
			int3[,5]=1-int3[,5]
			colnames(int3)=c("gene","T1","b1","t1","p")
			if(nrow(int3)!=0){	
				int3=data.frame(int3,e1="focus")
			}

			INT=rbind(int1,int2)
			INT=rbind(INT,int3)
			INT[,c(3,5)]=apply(INT[,c(3,5)],2,as.numeric)
			INT=INT[which(INT[,1]==gene),]
			if(min(na.omit(sign(INT$b1)))==max(na.omit(sign(INT$b1)))&max(na.omit(sign(INT$b1)))==(-s)){
				INT=format_int(INT)
				colnames(INT)=c("gene","T2","b2","t2","e2")
				INT=data.frame(gene=gene,T1=T1,b1=as.numeric(int[x,3]),t1=int[x,"t1"],e1=int[x,"e1"],INT[,2:5,drop=F])
				intout=rbind(intout,INT)
			}
		}
		out=rbind(out,intout)
	}
	return(out)	
}
