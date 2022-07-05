library(ggplot2)
draw=lava[which(!lava[,1] %in% red & !lava[,2] %in% red),]
head(draw)

loci$color=ifelse((loci[,2] %%2 )==0,"b","a")
draw$color=loci[draw$id,"color"]
draw$z=-qnorm(draw$p)*draw$rg/abs(draw$rg)

p=ggplot(draw,aes(id,z,color=color))+geom_point(alpha=0.5)+theme_minimal()+
scale_color_manual(breaks=c("a","b"),values=c("#003C67FF","#A73030FF"),guide=FALSE)
p=p+theme(axis.text.x=element_blank(),axis.title.x=element_blank())+ylab("local rg Z score")
p=p+geom_hline(yintercept=-2.64)+ylim(c(-16,20))
p=p+annotate(geom="text",label="Cholelithiasis-High cholesterol",x=291,y=-15,size=4)
p=p+annotate(geom="text",label="Atrial bibrillation-stroke",x=698,y=-11,size=4)
p=p+annotate(geom="text",label="Mouth Ulcers-Crohn's Disease",x=2100,y=-11,size=4)
p
p_f2a=p

pdf("f2a.pdf",height=3.5,width=10.5)
p_f2a
dev.off()

tiff("f2a.tiff",height=3,width=9,res=300,units="in")
p_f2a
dev.off()


####################f2b as corrplot###########################################
library(dplyr)
library(corrplot)
library(RColorBrewer)
d_f2a=loci_antag[which(!loci_antag[,1] %in% red & !loci_antag[,2] %in% red),c("pheno1","pheno2","id")]
int=data.frame(pheno1=d_f2a$pheno2,pheno2=d_f2a$pheno1,id=d_f2a$id)
d_f2a=rbind(d_f2a,int)

d_f2a=d_f2a %>% count(pheno1,pheno2)
d_f2a=data.frame(d_f2a)

d_f2a=spread(d_f2a,pheno1,n)
rownames(d_f2a)=d_f2a[,1]
d_f2a=d_f2a[,-1]
dim(d_f2a)
d_f2a[is.na(d_f2a)]=0
d_f2a[d_f2a<5]=0
d_f2a=d_f2a[ud,ud]
ud=which(rowSums(d_f2a)!=0)
,order="AOE"
rownames(d_f2a)=short[match(rownames(d_f2a),short[,1]),2]
colnames(d_f2a)=rownames(d_f2a)
 col = colorRampPalette(c("white", "#BEBADA"))(23)

pdf("figure 2b.pdf",height=6,width=6)
corrplot(as.matrix(d_f2a), is.corr = FALSE, type="lower",
tl.col = 'black', tl.srt = 45,tl.cex=0.8,method="circle")
dev.off()

g <- grab_grob()
 ptest=as.ggplot(g)

##############################################
l=table(loci_antag[which(!loci_antag[,1] %in% red & !loci_antag[,2] %in% red),"id"])
max(l)
which(l==25)

int=lava[which(lava$id==1353),]
ad=union(loci_antag[which(loci_antag$id==1353),1],loci_antag[which(loci_antag$id==1353),2])
setdiff(ad,red)

int=int[which(!int[,1] %in% red & !int[,2] %in% red),]
int$show=ifelse(int$fdr<0.05 & abs(int$rg)>0.5,1,0)
table(c(int[which(int$show==1),1],int[which(int$show==1),2]))
write.csv(int[which(int$show==1),],"edge.csv")

################################################
res=res[-which(res[,1] %in% red),c(1,3,4)]
dim(res)
summary(res[,2])
res=res[which(res[,2]>0.005),]

res$N=NA
for (i in 1:nrow(res)){
	dis=res[i,1]
	n=length(which(loci_antag[,1]==dis | loci_antag[,2]==dis))
	res[i,"N"]=n
}

res$name=abr[match(res[,1],abr[,1]),2]

library(ggplot2)
library(RColorBrewer)
library(gg.gap)
res=res[order(res[,3],decreasing=T),]
res$name=factor(res$name,levels=res$name[nrow(res):1])
colnames(res)=c("dis","propsnp","proph2","N.Tloci","name")
pf2d=ggplot(res,aes(name,proph2,size=N.Tloci))+geom_point(alpha=0.7,color=brewer.pal(3,"Set3")[3])+
theme_classic()+ylim(c(0,0.42))+scale_y_continuous(breaks=c(0,0.2,0.4),labels=c(0,20,60))+
theme(axis.ticks=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=11,
color="black"),axis.title.x=element_blank(),axis.title.y=element_text(color="black"),
legend.position=c(0.1,0.6))+labs(y="%h2 explained by Tloci",size="Number of Tloci") 

pdf("f2d.pdf",height=3,width=8.5)
pf2d
dev.off()

res=res[order(res[,6],decreasing=T),]
head(res)
head(res[which(res[,3]>0.005),])


gene_stat=data.frame(id=gene_stat[,1],n=as.numeric(gene_stat[,2]),mhc=gene_stat[,4],
gene=annot[match(gene_stat[,1],annot[,1]),2])

gene_stat=gene_stat[order(gene_stat[,2],decreasing=TRUE),]

head(gene_stat,n=20)

head(gene_stat[which(gene_stat[,3]==0),],n=50)

library(ggsci)
library(ggplot2)
draw$mhc=ifelse(draw$mhc=="yes","MHC","nonMHC")
draw=draw[order(draw$n,decreasing=TRUE),]
draw$gene=factor(draw$gene,levels=draw[37:1,3])
p=ggplot(draw,aes(gene,n,color=mhc,fill=mhc))+geom_bar(stat="identity",width=1)+
facet_grid(mhc~.,scales="free_y",space="free")+coord_flip()+theme_classic()
p=p+theme(legend.position="non",axis.ticks=element_blank(),axis.text=element_text(color="black"),
axis.line.y=element_blank())+ylab("N.disease associated")+scale_y_continuous(expand = c(0,0))+
scale_color_manual(breaks=c("MHC","nonMHC"),values=c("#003C67FF","#A73030FF"),guide=FALSE)+
scale_fill_manual(breaks=c("MHC","nonMHC"),values=c("#003C67FF","#A73030FF"),guide=FALSE)
p_f4b=p
save(p_f4b,file="p_f4b.RData")
pdf("f4b.pdf",height=7,width=3.5)
p_f4b
dev.off()

disname=read.csv("disname.csv",head=F)
d_f4a=data.frame(from=disname[match(dis_network[,1],disname[,1]),2],to=
disname[match(dis_network[,2],disname[,1]),2],n=as.numeric(dis_network[,3]))
d_f4a=d_f4a[which(!is.na(d_f4a[,1]) & !is.na(d_f4a[,2]) & d_f4a[,3]!=0),]
d_f4a=d_f4a[which(d_f4a$n>100),]
d_f4a=spread(d_f4a,from,n)
d_f4a[is.na(d_f4a)]=0
rownames(d_f4a)=d_f4a[,1]
d_f4a=d_f4a[,-1]
circos.clear()

pdf("p_f4a.pdf",height=5,width=5)
ord=c("CD","Hypothyroidism","Asthma","SCZ","glaucoma","HBP","skin cancer","CAD","High cholesterol",
"UC","Vitiligo","T2D","BCC","eczema","GERD","myopia",
"gout","ulcer","RA","JIA","PSC","PBC","AF")

col_fun = colorRamp2(c(50,329), c("#FFEEEE", "#FF0000"), transparency = 0.5)
chordDiagram(as.matrix(d_f4a),order = ord, annotationTrack = c("name","grid"),
annotationTrackHeight=c(0.03,0.07))
dev.off()

write.csv(d_f4a,"edge.csv",quote=F)

#############################f4a as corr plot###############################
dis_network=dis_network[which(dis_network[,1] %in% short[,1]),]
d_f4a=data.frame(from=short[match(dis_network[,1],short[,1]),2],to=
short[match(dis_network[,2],short[,1]),2],n=as.numeric(dis_network[,4]))
d_f4a=d_f4a[which(!is.na(d_f4a[,1]) & !is.na(d_f4a[,2]) & d_f4a[,3]!=0),]
d_f4a=d_f4a[which(d_f4a$n>20),]
int=data.frame(from=d_f4a$to,to=d_f4a$from,n=d_f4a$n)
d_f4a=rbind(d_f4a,int)
d_f4a=spread(d_f4a,from,n)
rownames(d_f4a)=d_f4a[,1]
d_f4a=d_f4a[,-1]
dim(d_f4a)
d_f4a[is.na(d_f4a)]=0

pdf("f4a.pdf",height=7,width=7)
 col = colorRampPalette(c("#67001F","#FFFFFF", "#053061"))(201)
corrplot(as.matrix(d_f4a), is.corr = FALSE, type="lower",
tl.col = 'black', tl.srt = 45,tl.cex=0.8,method="circle",order="AOE")
dev.off()


########################f4c heatmap######################
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
gene="ENSG00000244731"

dis=all_assoc[which(all_assoc[,1]==gene),2]

res_smr=all_smr[which(all_smr[,1]==gene & all_smr$trait %in% dis),]
res_smr$z=qnorm(res_smr[,4])*as.numeric(res_smr[,3])/abs(as.numeric(res_smr[,3]))
res_smr=res_smr[order(res_smr[,4]),]
res_smr=res_smr[!duplicated(res_smr$trait),c("trait","z")]

res_focus=focusres[which(focusres[,1]==gene & focusres$trait %in% dis),]
res_focus=res_focus[order(res_focus[,5]),]
res_focus=res_focus[!duplicated(res_focus$trait),c("trait","twas_z")]

res_pred=allpred[which(allpred[,1]==gene & allpred$trait %in% dis),]
res_pred$z=qnorm(res_pred[,4])*as.numeric(res_pred[,3])/abs(as.numeric(res_pred[,3]))
res_pred=res_pred[order(res_pred[,4]),]
res_pred=res_pred[!duplicated(res_pred$trait),c("trait","z")]

d_f4c=data.frame(smr=res_smr[match(dis,res_smr[,1]),2],pred=res_pred[match(dis,res_pred[,1]),2],
focus=res_focus[match(dis,res_focus[,1]),2])
rownames(d_f4c)=dis
d_f4c[is.na(d_f4c)]=0
d_f4c[,1:2]=-d_f4c[,1:2]
x=disname[match(rownames(d_f4c),disname[,1]),2]
d_f4c=d_f4c[!is.na(x),]
rownames(d_f4c)=na.omit(x)
colnames(d_f4c)=c("SMR","SPredixcan","FOCUS")
col_fun= colorRamp2(c(-20,0,20),brewer.pal(11,"RdBu")[c(11,6,1)])

pdf("f4c.pdf",height=2,width=6.2)
Heatmap(t(df4c),col=col_fun,show_row_dend=F,show_column_dend=F,heatmap_legend_param = 
list(at = c(-20,0,20),title="C4A Z score",direction = "horizontal",
border="black",title_position="topcenter"),
row_names_gp = gpar(fontsize =10),column_names_gp = gpar(fontsize =10))
dev.off()

gene="ENSG00000108379"

dis=all_assoc[which(all_assoc[,1]==gene),2]

res_smr=all_smr[which(all_smr[,1]==gene & all_smr$trait %in% dis),]
res_smr$z=qnorm(res_smr[,4])*as.numeric(res_smr[,3])/abs(as.numeric(res_smr[,3]))
res_smr=res_smr[order(res_smr[,4]),]
res_smr=res_smr[!duplicated(res_smr$trait),c("trait","z")]

res_focus=focusres[which(focusres[,1]==gene & focusres$trait %in% dis),]
res_focus=res_focus[order(res_focus[,5]),]
res_focus=res_focus[!duplicated(res_focus$trait),c("trait","twas_z")]

res_pred=allpred[which(allpred[,1]==gene & allpred$trait %in% dis),]
res_pred$z=qnorm(res_pred[,4])*as.numeric(res_pred[,3])/abs(as.numeric(res_pred[,3]))
res_pred=res_pred[order(res_pred[,4]),]
res_pred=res_pred[!duplicated(res_pred$trait),c("trait","z")]

d_f4c=data.frame(smr=res_smr[match(dis,res_smr[,1]),2],pred=res_pred[match(dis,res_pred[,1]),2],
focus=res_focus[match(dis,res_focus[,1]),2])
rownames(d_f4c)=dis
d_f4c[is.na(d_f4c)]=0
d_f4c[,1:2]=-d_f4c[,1:2]
x=disname[match(rownames(d_f4c),disname[,1]),2]
d_f4c=d_f4c[!is.na(x),]
rownames(d_f4c)=na.omit(x)
colnames(d_f4c)=c("SMR","SPredixcan","FOCUS")
col_fun= colorRamp2(c(-10,0,10),brewer.pal(11,"RdBu")[c(11,6,1)])


ph=Heatmap(t(d_f4c),col=col_fun,show_row_dend=F,show_column_dend=F,heatmap_legend_param = 
list(at = c(-10,0,10),title="WNT3 Z score",legend_title_gp=gpar(fontsize=6),
border="black",direction = "horizontal",title_position="topcenter"),
row_names_gp = gpar(fontsize =10),column_names_gp = gpar(fontsize =10),cluster_rows=F)
pdf("f4d.pdf",height=2.4,width=5)
draw(ph,heatmap_legend_side="bottom")
dev.off()

######################word cloud##############################
library(clusterProfiler)
library(org.Hs.eg.db)

Agene=unique(all_assoc[,1])
Tgene=unique(antagres[,1])
Ngene=setdiff(Agene,Tgene)

gene.d =bitr(Ngene, fromType = "ENSEMBL", 
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db)
gl=as.character(gene.d[,2])
NGO=enrichGO(gene         =gl,
                OrgDb         = org.Hs.eg.db,
                
                ont           = "BP",
                pAdjustMethod = "BH",
              
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
NGO=data.frame(simplify(NGO))
head(TGO)
NGO[1:30,c(2,6)]
write.csv(NGO,"NGO.csv")

source("https://gist.githubusercontent.com/jokergoo/bfb115200df256eeacb7af302d4e508e/raw/14f315c7418f3458d932ad749850fd515dec413b/word_cloud_grob.R")
getR=function(x){
x1=unlist(strsplit(as.character(x[1]),"/",fixed=TRUE))
x1=as.numeric(x1[1])/as.numeric(x1[2])
x2=unlist(strsplit(as.character(x[2]),"/",fixed=TRUE))
x2=as.numeric(x2[1])/as.numeric(x2[2])
r=x1/x2
return(r)
}

library(stringr)
words=TGO[c("GO:0060333","GO:0050852","GO:0019058","GO:0006914","GO:0002262","GO:0051260",
"GO:0048002"),c(2)]
pval=-log10(TGO[c("GO:0060333","GO:0050852","GO:0019058","GO:0006914","GO:0002262","GO:0051260",
"GO:0048002"),c(6)])

r=apply(TGO[c("GO:0060333","GO:0050852","GO:0019058","GO:0006914","GO:0002262","GO:0051260",
"GO:0048002"),3:4],1,getR)
words=str_wrap(words, width = 30)


col_fun= colorRamp2(c(1.5,5),brewer.pal(11,"RdBu")[c(11,1)])

gb_i = word_cloud_grob(words, fontsize = 3*r, max_width = unit(30, "mm"),
col=col_fun(pval))

pdf("f4e.pdf",width=6,height=5)
grid.newpage()
grid.draw(gb_i)
dev.off()

########### f4f ###############

words=NGO[c("GO:0006644","GO:0006066","GO:0009150","GO:0043087","GO:0044782","GO:0007265","GO:0032147"),c(2)]
pval=-log10(NGO[c("GO:0006644","GO:0006066","GO:0009150","GO:0043087","GO:0044782","GO:0007265","GO:0032147"),c(6)])

r=apply(NGO[c("GO:0006644","GO:0006066","GO:0009150","GO:0043087","GO:0044782","GO:0007265","GO:0032147"),3:4],1,getR)
words=str_wrap(words, width = 37)


col_fun= colorRamp2(c(1.5,7),brewer.pal(11,"RdBu")[c(11,1)])

gb_i = word_cloud_grob(words, fontsize = 3*r, max_width = unit(20, "mm"),
col=col_fun(pval))

pdf("f4f.pdf",width=6,height=6)
grid.newpage()
grid.draw(gb_i)
dev.off()

########### f4g ###############
head(antagres)
dp=c("ID503")
hgene=antagres[which(antagres$T1 == dp | antagres$T2 == dp),1]
setdiff(hgene,mhcgene)
gene.d =bitr(hgene, fromType = "ENSEMBL", 
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db)
gl=as.character(gene.d[,2])
hGO=enrichGO(gene         =gl,
                OrgDb         = org.Hs.eg.db,
                
                ont           = "BP",
                pAdjustMethod = "BH",
              
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
hGO=data.frame(simplify(hGO))
head(hGO)
NGO[1:30,c(2,6)]
write.csv(HGO,"HGO.csv")

words=HGO[c("GO:1903749","GO:0019882","GO:1901660"),]
pval=-log10(eGO[c("GO:1903749","GO:0019882","GO:1901660"),c(6)])

r=apply(eGO[c("GO:1903749","GO:0019882","GO:1901660"),3:4],1,getR)
words=str_wrap(words, width = 40)


col_fun= colorRamp2(c(1.3,1.6),brewer.pal(11,"RdBu")[c(1,11)])

gb_i = word_cloud_grob(words, fontsize = 3*r, max_width = unit(50, "mm"),
col=col_fun(pval))

pdf("f4g.pdf",width=6,height=6)
grid.newpage()
grid.draw(gb_i)
dev.off()

########### f4h ###############
head(antagres)
dp=c("ID285")
egene=antagres[which(antagres$T1 == dp | antagres$T2 == dp),1]

gene.d =bitr(egene, fromType = "ENSEMBL", 
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db)
gl=as.character(gene.d[,2])
eGO=enrichGO(gene         =gl,
                OrgDb         = org.Hs.eg.db,
                
                ont           = "BP",
                pAdjustMethod = "BH",
              
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
eGO=data.frame(simplify(eGO))
head(eGO)
NGO[1:30,c(2,6)]

words=eGO[c("GO:0002495","GO:0007159","GO:0050870"),c(2)]
pval=-log10(eGO[c("GO:0002495","GO:0007159","GO:0050870"),c(6)])

r=apply(eGO[c("GO:0002495","GO:0007159","GO:0050870"),3:4],1,getR)
words=str_wrap(words, width = 37)


col_fun= colorRamp2(c(2,3.2),brewer.pal(11,"RdBu")[c(11,1)])

gb_i = word_cloud_grob(words, fontsize = 3*r, max_width = unit(50, "mm"),
col=col_fun(pval))

pdf("f4h.pdf",width=6,height=6)
grid.newpage()
grid.draw(gb_i)
dev.off()
