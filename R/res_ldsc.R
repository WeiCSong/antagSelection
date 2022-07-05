res_gene_all$h2=h2[res_gene_all[,5],2]
res_gene_all=res_gene_all[!is.na(res_gene_all$h2),]
res_gene_all$n=n[res_gene_all[,5],4]
res_gene_all=res_gene_all[!is.na(res_gene_all$n),]

d=res_gene_all[which(!res_gene_all$disease %in% red & 
res_gene_all[,1]=="antaggene" & res_gene_all$h2<0.2),]
mgen=metagen(d[,2]*1e7/d[,6],d[,3]*1e7/d[,6])
mgen=update.meta(mgen, prediction = TRUE)
mgen

#qtl_bg_all:1.76,18.87

d$tao=d[,2]*d[,7]/d[,6]
plot(d[,2],d$h2)

ivw=function(x,y){
	b=sum(x/y)/sum(1/y)
	se=1/sum(1/y)
	return(c(b,se))
}

ivw(d[,2]*d[,7]/d[,6],d[,3]*d[,7]/d[,6])

###gene result
#all vs bg, h2>0.2: 1.071+-0.015
#all vs bg, h2<0.2: 1.271+-0.004
#ple vs bg, h2>0.2: 1.151+-0.016
#ple vs bg, h2<0.2: 1.247+-0.005
#ant vs bg, h2>0.2: 1.212+-0.020
#ant vs bg, h2<0.2: 1.058+-0.005
#ple vs all,h2>0.2: 0.228+-0.035
#ple vs all,h2<0.2: 0.226+-0.012
#ant vs all,h2>0.2: 0.072+-0.032
#ant vs all,h2<0.2: -0.112+-0.010

d=res_qtl_all[which(!res_qtl_all$disease %in% red & 
res_qtl_all[,1]=="antaggene" & res_qtl_all$h2>0.2),]
ivw(d[,2]*d[,7]/d[,6],d[,3]*d[,7]/d[,6])

###qtl result
ldsc=read.csv("res_ldsc_loci.csv")
library(ggplot2)
d=ldsc[which(ldsc$type=="tva"),]




ldsc$h2=h2[ldsc[,2],2]
ldsc=ldsc[!is.na(ldsc$h2),]
ldsc$n=n[ldsc[,2],4]
ldsc=ldsc[!is.na(ldsc$n),]
colnames(ldsc)[1:4]=c("disease","type","b","z")
ldsc$se=ldsc$b/ldsc$z

res_qtl_all$tao=res_qtl_all[,2]*res_qtl_all[,7]/res_qtl_all[,6]
res_qtl_all$taose=res_qtl_all[,3]*res_qtl_all[,7]/res_qtl_all[,6]
res_qtl_bg$tao=res_qtl_bg[,2]*res_qtl_bg[,7]/res_qtl_bg[,6]
res_qtl_bg$taose=res_qtl_bg[,3]*res_qtl_bg[,7]/res_qtl_bg[,6]

res_qtl_all$disname=disname[res_qtl_all$disease,2]
res_gene_all$disname=disname[res_gene_all$disease,2]
res_qtl_bg$disname=disname[res_qtl_bg$disease,2]
res_gene_bg$disname=disname[res_gene_bg$disease,2]

write.csv(res_qtl_all,"qtl_ldsc.csv")
write.csv(res_gene_all,"gene_ldsc.csv")
write.csv(res_qtl_bg,"bqtl_ldsc.csv")
write.csv(res_gene_bg,"bgene_ldsc.csv")

head(res_gene_bg)
int=spread(res_gene_bg[,c("Name","tao","disease")],Name,tao)
t.test(int[,2],int[,3],paired=TRUE)

head(res_gene_all)

############figure 3###########################
d=ldsc[which(!ldsc$disease %in% red & ldsc$type=="tva" & ldsc$h2>0.2),]
d$tao=d$b*d$n/d$h2
d$taose=d$se*d$n/d$h2
d$dis=disname[match(d$disease,disname[,1]),2]

ivw=function(x,y){
	b=sum(x/y)/sum(1/y)
	se=1/sum(1/y)
	return(c(b,se))
}

ivw(d$tao,d$taose)
draw=rbind(d[,8:10],data.frame(tao=0.12,taose=0.01,dis="meta"))
draw$dis=factor(draw$dis,levels=draw$dis[c(10,9,3,2,8,1,6,4,5,7)])
labels=levels(draw$dis)
labels[3]="ADHD"
labels[2]="thyroid cancer"
labels[8]="ASD"
labels[10]="PSC"

p=ggplot(draw,aes(tao,dis,color=dis))+geom_point(size=2)+theme_minimal()+
geom_errorbarh(aes(xmax=tao+1.96*taose,xmin=tao-1.96*taose),height=0)+
scale_color_manual(breaks=levels(draw$dis),values=ifelse(levels(draw$dis)==
"meta","red","black"))+theme(legend.position="none",axis.title.y=element_blank())+
geom_vline(xintercept=0)+scale_y_discrete(labels= labels)+xlab("Tloci vs Aloci 而")
forest_loci_ldsc=p

d=ldsc[which(!ldsc$disease %in% red & !ldsc$type=="tva"),]
d$tao=d$b*d$n/d$h2
d$taose=d$se*d$n/d$h2
d$dis=disname[match(d$disease,disname[,1]),2]

summary(d[which(d[,2]=="all"),"tao"])

int=spread(d[,c("disease","type","z")],type,z)
t.test(int[,3],int[,2],paired=T)


d=split(d,f=d$type)
draw=d[[2]][,8:10]
colnames(draw)=c("tt","tse","dis")
draw$at=d[[1]][match(draw$dis,d[[1]]$dis),"tao"]
draw$ase=d[[1]][match(draw$dis,d[[1]]$dis),"taose"]
draw=draw[-which(draw$dis %in% disname[match(red,disname[,1]),2]),]
p=ggplot(draw,aes(at,tt))+geom_point(size=1)+theme_minimal()+geom_errorbar(
aes(ymax=tt+1.96*tse,ymin=tt-1.96*tse),size=0.3,width=0)+geom_errorbarh(aes(xmax=
at+1.96*ase,xmin=at-1.96*ase),size=0.3,height=0)+geom_segment(x=0,y=0,xend=6,yend=6)+
xlab("Aloci 而")+ylab("Tloci 而")
dot_ldsc_loci=p

p_f3b=ggplot(d,aes(x=z))+geom_density(fill="#003C67",alpha=0.5)+theme_classic()+
geom_vline(xintercept=0)+xlab("Z of Tloci vs Aloci 而")

###############figure 3d####################
draw=bfull[which(bfull[,1]=="B2"),]
which(draw[,2]>0)

p_f3d=ggplot(draw,aes(x=X2))+geom_density(fill="#003C67",alpha=0.5)+theme_classic()+
geom_vline(xintercept=0)+xlab("t of Tloci vs Nloci")+annotate(geom="text",label=
"B2 statistics\ncomparison",size=2.5,x=3.2,y=0.2)



t.test(bfull[which(bfull[,1]=="B2"),2])

draw=bfull[which(bfull[,1]=="B2" & bfull[,3] %in% hh),]


draw[,3]=disname[match(draw[,3],disname[,1]),2]
draw=draw[,-1]
draw=rbind(draw,data.frame(X2=mean(draw[,1]),dis="meta"))
draw[3,2]="ASD"
draw[4,2]="ADHD"
draw[7,2]="PSC"

colnames(draw)=c("t","disease")
draw$disease=factor(draw$disease,levels=draw[9:1,2])

p=ggplot(draw,aes(t,disease,color=disease))+geom_point(size=2)+theme_minimal()+
geom_errorbarh(aes(xmax=t+1.96*se,xmin=t-1.96*se),height=0)+
scale_color_manual(breaks=levels(draw$disease),values=ifelse(levels(draw$disease)==
"meta","red","black"))+theme(legend.position="none")+ylab("")+
geom_vline(xintercept=0)+xlab("t of Tloci vs Nloci")
forest_loci_b2=p

draw$se=c(1,1,1,1,1,1,1,1,0.125)

x=int1[is.finite(int1$B2),]
max(x[which(x$group=="no"),"B2"])


################################figure 3f#########################
library(tidyr)
library(dplyr)
head(prs)
prs$group=paste(prs$study,prs$type,sep="_")
prs$group=paste(prs$group,prs$source,sep="_")
prsd=spread(prs[,c(1,6,7)],group,p)
draw=na.omit(prs[which(prs$study=="loci"&prs$source=="relate"),c(1,3,6)])
draw[,3]=ifelse(draw[,3]<0.001,"yes","no")
draw[,2]=ifelse(draw[,2]=="non","Nloci","Tloci")

fisher.test(cbind(table(draw[which(draw[,2]=="Nloci"),3]),table(draw[which(draw[,2]=="Tloci"),3])))
colnames(draw)[3]=c("p")


p_f3f=ggplot(draw,aes(type,fill=p,color=p))+geom_bar()+theme_minimal()+
theme(axis.title.y=element_blank(),axis.text.x=element_text(color="black",size=10),
legend.margin=margin(t=0, r=-0.6, b=0, l=-0.2, unit="cm"))+
scale_fill_manual(breaks=c("yes","no"),values=c("#A73030FF","#003C67FF"),name=
"Fisher p\n=0.006\n\nsignificant\nadaptation")+xlab("PRS group")+
scale_color_manual(breaks=c("yes","no"),values=c("#A73030FF","#003C67FF"),guide=F)

####################traj##########################
p_scz=p_scz+geom_point(size=2)+theme(legend.position=c(0.3,0.2))
p_jia=p_f3g+geom_point(size=2)+theme(legend.position=c(0.3,0.2))
p_eczema=p_eczema+geom_point(size=2)+theme(legend.position=c(0.3,0.2))
p_prostate=p_prostate+geom_point(size=2)+theme(legend.position=c(0.3,0.2))

f3up=dot_ldsc_loci+p_f3b+forest_loci_ldsc+plot_layout(width=c(3,2,3))
f3mid=p_f3d+forest_loci_b2
f3down=p_f3f+guide_area()+p_scz+p_prostate+plot_layout(width=c(1.5,0.1,3.2,3.2))
f3=f3up/f3mid/f3down+plot_annotation(tag_levels='A')

CairoPDF("Figure 3.pdf",height=10,width=8)
f3
dev.off()

tiff("Figure traj.tiff",res=300,units="in",height=3,width=6)
p_eczema+p_jia
dev.off()


##################figure 5##############
res_gene_all$tao=res_gene_all[,2]*res_gene_all[,7]/res_gene_all[,6]
res_gene_all$taose=res_gene_all[,3]*res_gene_all[,7]/res_gene_all[,6]

res_gene_bg$tao=res_gene_bg[,2]*res_gene_bg[,7]/res_gene_bg[,6]
res_gene_bg$taose=res_gene_bg[,3]*res_gene_bg[,7]/res_gene_bg[,6]

d1=spread(res_gene_bg[,c("Name","disease","tao")],Name,tao)
d2=spread(res_gene_bg[,c("Name","disease","taose")],Name,taose)
draw=data.frame(dis=d1$disease,all=d1$allgene,allse=d2$allgene,antag=d1$antaggene,
antagse=d2$antaggene)

library(ggplot2)
p_f5a=ggplot(draw,aes(all,antag))+geom_point()+geom_errorbar(aes(ymax=antag+1.96*antagse,
ymin=antag-1.96*se),width=0)+geom_errorbarh(aes(xmax=all+1.96*allse,xmin=all-1.96*allse),
height=0)+theme_minimal()+xlim(c(-2,10))+ylim(c(-1,10))+xlab("Agene \u03c4")+ylab("Tgene \u03c4")+
geom_segment(x=0,y=0,xend=10,yend=10)

draw=rbind(res_gene_all[which(res_gene_all[,5] %in% hh & res_gene_all[,1]=="antaggene"),c(8,9,5)],data.frame(disease="meta",tao=0.07,taose=0.03))
draw$disease=labels[10:1]
draw$disease=factor(draw$disease,levels=labels)
labels=levels(draw$dis)

labels[3]="glaucoma"
labels[2]="thyroid cancer"
labels[8]="ASD"
labels[4]="PSC"
labels[5]="Schizophrenia"
labels[6]="Ulcerative colitis"
labels[7]="Bipolar disorder"
labels[9]="Crohn's disease"
labels[10]="ADHD"


p=ggplot(draw,aes(tao,dis,color=dis))+geom_point(size=3)+theme_minimal()+
geom_errorbarh(aes(xmax=tao+1.96*taose,xmin=tao-1.96*taose),height=0)+
scale_color_manual(breaks=levels(draw$dis),values=ifelse(levels(draw$dis)==
"meta","red","black"))+theme(legend.position="none")+ylab("")+
geom_vline(xintercept=0)+scale_y_discrete(labels= labels)+xlab("Tgene vs Agene \u03c4")
forest_gene_ldsc=p

d=res_gene_all[which(res_gene_all[,1]=="antaggene"),c(8,9,5)]
d$z=d[,1]/d[,2]
p_f5b=ggplot(d,aes(x=z))+geom_density(fill="#003C67",alpha=0.5)+theme_classic()+
geom_vline(xintercept=0)+xlab("Z of Tgene vs Agene \u03c4")
#######################figure 5 b2###########################
draw=genecm[which(genecm[,1]=="B2"),]
p_f5d=ggplot(draw,aes(x=X2))+geom_density(fill="#003C67",alpha=0.5)+theme_classic()+
geom_vline(xintercept=0)+xlab("t of Tgene vs Ngene")+annotate(geom="text",label=
"B2 statistics\ncomparison",size=2.5,x=1.5,y=0.45)

draw=genecm[which(genecm[,1]=="B2" & genecm[,3] %in% hh),]
draw[,3]=disname[match(draw[,3],disname[,1]),2]
draw=draw[,-1]
draw=rbind(draw,data.frame(X2=mean(draw[,1]),dis="meta"))

draw[5,2]="ADHD"
draw[7,2]="PSC"

colnames(draw)=c("t","disease")
draw=draw[c(7:1,8),]
draw$disease=factor(draw$disease,levels=draw[8:1,2])
draw$se=c(1,1,1,1,1,1,1,1/7)
p=ggplot(draw,aes(t,disease,color=disease))+geom_point(size=3)+theme_minimal()+
geom_errorbarh(aes(xmax=t+1.96*se,xmin=t-1.96*se),height=0)+
scale_color_manual(breaks=levels(draw$disease),values=ifelse(levels(draw$disease)==
"meta","red","black"))+theme(legend.position="none",axis.title.y=element_blank())+
geom_vline(xintercept=0)+xlab("t of Tgene vs Ngene")
forest_gene_b2=p

draw=na.omit(prs[which(prs$study=="gene"&prs$source=="relate"),c(1,3,6)])
draw[,3]=ifelse(draw[,3]<0.001,"yes","no")
draw[,2]=ifelse(draw[,2]=="non","Nloci","Tloci")

fisher.test(cbind(table(draw[which(draw[,2]=="Nloci"),3]),table(draw[which(draw[,2]=="Tloci"),3])))

colnames(draw)[3]=c("p")


p_f5f=ggplot(draw,aes(type,fill=p,color=p))+geom_bar()+theme_minimal()+
theme(axis.title.y=element_blank(),axis.text.x=element_text(color="black",size=10),
legend.margin=margin(t=0, r=-0.5, b=0, l=0, unit="cm"))+scale_fill_manual(
breaks=c("yes","no"),values=c("#A73030FF","#003C67FF"),name=
"Fisher p\n=0.62\n\nsignificant\nadaptation")+xlab("PRS group")+
scale_color_manual(breaks=c("yes","no"),values=c("#A73030FF","#003C67FF"),guide=F)

p_f5f=p_f5f+scale_x_discrete(labels=c("Ngene","Tgene"))
library(patchwork)
f5up=p_f5a+p_f5b+forest_gene_ldsc+plot_layout(width=c(3,2,3))
f5down=p_f5d+forest_gene_b2+p_f5f+guide_area()+plot_layout(width=c(3,3,1.9,0.1))
f5=f5up/f5down+ plot_annotation(tag_levels = 'A')

CairoPDF("Figure 5.pdf",height=7,width=8)
f5
dev.off()

