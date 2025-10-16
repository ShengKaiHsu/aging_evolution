rm(list=ls())
library(ExpressionNormalizationWorkflow)
library(edgeR)
library(pheatmap)
library(biomaRt)
library(pathview)
library(pvca)
library(sva)
library(topGO)
library(RcisTarget)
library(qqman)
library(VennDiagram)
library(scales)
library(maSigPro)
library(WGCNA)
setwd("/Volumes/Temp2/shengkai/age_CGE/formal_experiment/")
####import basic function####
cont_table=function(query,background,classifyer){
  p1=length(Reduce(intersect,list(query,background,classifyer)))
  q1=length(intersect(query,background))-p1
  p0=length(setdiff(intersect(background,classifyer),intersect(query,classifyer)))
  q0=length(setdiff(background,query))-p0
  return(matrix(c(p1,p0,q1,q0),2,2))
}
genomic_enrichment=function(query,background,ann_summary,windowsize,steps,chr){
  X0=0
  X=windowsize
  S=steps
  pos=c()
  counts=c()
  fold_enrichment=c()
  p.val=c()
  x=ann_summary[ann_summary[,1]%in%chr,]
  while (X0<max(x[,c(2,3)])){
    gene_within=x[x[,2]<X0+X/2&x[,2]>X0-X/2,5]
    p=c()
    fe=c()
    for (j in 1:length(query)){
      p=cbind(p,fisher.test(cont_table(query[[j]],background,gene_within),alternative = "greater")$p.value) 
      fe=cbind(fe,(length(intersect(gene_within,query[[j]]))/length(gene_within))/(length(query[[j]])/length(background)))
    }
    p.val=rbind(p.val,p)
    fold_enrichment=rbind(fold_enrichment,fe)
    colnames(p.val)=names(query)
    colnames(fold_enrichment)=names(query)
    counts=c(counts,length(x[x[,2]<X0+X/2&x[,2]>X0-X/2,5]))
    pos=c(pos,X0)
    X0=X0+S
  }
  out_p=data.frame("mid_pos"=pos,"total_gene_counts"=counts,p.val)
  out_fe=data.frame("mid_pos"=pos,"total_gene_counts"=counts,fold_enrichment)
  out=list("p"=out_p,"fe"=out_fe)
  return(out)
}
ID_converter=function(ID,db,attributes,filters){
  getBM(attributes=attributes,filters=filters,mart=db,values=ID)
}
ROC_like_curve=function(a,b,bin,mycol="black",positive=T){
  lvl=seq(0,1,bin)
  if (positive==T){
    inter=c()
    for (i in lvl){
      conserved_genes=intersect(rownames(a),rownames(b))
      tmp1=a[conserved_genes,]
      tmp2=b[conserved_genes,]
      index=sum(rank(tmp1$logFC,ties.method = "random")>length(conserved_genes)*(1-i)&rank(tmp2$logFC,ties.method = "random")>length(conserved_genes)*(1-i))
      inter=c(inter,index/(length(conserved_genes)*i))
    }
    print(data.frame(r=cor(tmp1$logFC,tmp2$logFC),AUC=auc(lvl,inter,0,1,type = "spline")))
    points(lvl,inter,asp=1,type="l",lwd=3,col=mycol)
  }
  else {
    inter=c()
    for (i in lvl){
      conserved_genes=intersect(rownames(a),rownames(b))
      tmp1=a[conserved_genes,]
      tmp2=b[conserved_genes,]
      index=sum(rank(tmp1$logFC,ties.method = "random")<length(conserved_genes)*i&rank(tmp2$logFC,ties.method = "random")<length(conserved_genes)*i)
      inter=c(inter,index/(length(conserved_genes)*i))
    }
    print(data.frame(r=cor(tmp1$logFC,tmp2$logFC),AUC=auc(lvl,inter,0,1,type = "spline")))
    points(lvl,inter,asp=1,type="l",lwd=3,col=mycol)
  }
}
ROC_like_curve_GO=function(a,b,bin,mycol="black"){
  lvl=seq(0,1,bin)
  inter=c()
  for (i in lvl){
    conserved_GO=intersect(a$GO.ID,b$GO.ID)
    tmp1=a[a$GO.ID%in%conserved_GO,]
    tmp2=b[b$GO.ID%in%conserved_GO,]
    index=length(intersect(tmp1[0:(length(conserved_GO)*i),1],tmp2[0:(length(conserved_GO)*i),1]))
    inter=c(inter,index/(length(conserved_GO)*i))
  }
  print(auc(lvl,inter,0,1,type = "spline"))
  points(lvl,inter,asp=1,type="l",lwd=3,col=mycol)
}
pseudo_chr=function(x){
  for (i in 1:length(unique(x$CHR))){
    x$pseudo_CHR[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]=i
  }
  return(x)
}
pseudo_pos=function(x){
  for (i in 1:length(unique(x$CHR))){
    if (i==1) x$pseudo_POS[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]=x$BP[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]
    else x$pseudo_POS[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]=x$BP[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]+max(x$pseudo_POS[x$CHR%in%c("X","2L","2R","3L","3R","4")[i-1]])
  }
  return(x)
}
p_resolution=function(x){
  x$P=10^(-as.numeric(strsplit2(x$"-logp",split = "=")[,2]))
  return(x)
}
SNPID_gen=function(x){
  x$SNPID=paste(x$CHR,x$BP,sep = "_")
  return(x)
}
snp_identifier=function(query,snpset,threshold1,threshold2){
  subset1=snpset[snpset$CHR%in%query$Chr&snpset$CHR=="X"&snpset$P<threshold1,]
  subset2=snpset[snpset$CHR%in%query$Chr&snpset$CHR!="X"&snpset$P<threshold2,]
  subset=rbind(subset1,subset2)
  count=sum(query$start<subset$BP&query$end>subset$BP,na.rm = T)
  SNPID=as.character(na.omit(subset$SNPID[query$start<subset$BP&query$end>subset$BP]))
  
  query$SNP_counts=count
  if(length(SNPID)>0) query$SNP_ID=paste(SNPID,collapse = ",")
  else query$SNP_ID=NA
  return(query)
}

####phenotype: mortality rate####
suv_rate=read.csv("../surv_updated.csv",header = T,row.names = 1,stringsAsFactors = F)
mort_rate=1-suv_rate
mort_rate_sample=mort_rate[,seq(1,21,2)]
mort_rate_sample[,1]=0
colnames(mort_rate_sample)=substr(colnames(mort_rate_sample),2,3)
mort_rate_sample=melt(mort_rate_sample)
mort_rate_sample$evo=rep(c("B","H"),each=6)
mort_rate_sample$rep=sprintf("%02d", 1:12)
mort_rate_sample$group=paste(paste0(mort_rate_sample$evo,mort_rate_sample$rep),mort_rate_sample$variable,sep = "_")
rownames(mort_rate_sample)=mort_rate_sample$group  
mort_rate_sample$scale=log(mort_rate_sample$value+1)/max(log(mort_rate_sample$value+1),na.rm = T)
idx=seq(0,0.8,0.05)
mort_rate_sample$sudo_age=c()
for (i in 1:16){
  mort_rate_sample$sudo_age[mort_rate_sample$value>=idx[i]&mort_rate_sample$value<idx[i+1]]=i
}
mort_rate_sample$sudo_age[is.na(mort_rate_sample$sudo_age)]=c(1,1,1,2,2,3)
  
avg_suv_rate=apply(suv_rate,2,function(x) tapply(x,gl(2,6,12), function(y) mean(y,na.rm=T)))
avg_suv_rate_sample=avg_suv_rate[,seq(1,21,2)]
avg_suv_rate_sample[,1]=1
avg_mort_rate_sample=1-avg_suv_rate_sample

#png("survival_rate_updated.png",width = 8.7,height = 6,units = "cm",pointsize = 6,res=600)
par(mar=c(5,5,2,2))
boxplot(suv_rate[1:6,-1:-2], col=alpha("forestgreen",0.5),border="Forestgreen",
        at=c(1:8,seq(10,18,2),seq(23,58,5))-0.25, 
        ylim=c(0,1),xlim=c(0,60),xaxt="n",cex=0.4,boxwex=0.4,xlab="Age (day)",ylab="Survial rate")
boxplot(suv_rate[-1:-6,-1:-2], add=T,col=alpha("salmon",0.5),border="salmon",
        at=c(1:8,seq(10,18,2),seq(23,58,5))+0.25,axes=F,boxwex=0.4,cex=0.4)
axis(1,at=c(1:8,seq(10,18,2),seq(23,58,5)))
#dev.off()

####phenotype: mortality rate (including cold evolved)####
suv_rate2=read.csv("../surv_rate_new.csv",header = T,row.names = 1,stringsAsFactors = F)

#png("../survival_rate_newCGE.png",width = 8.7,height = 6,units = "cm",pointsize = 6,res=600)
par(mar=c(5,5,2,2))
boxplot(suv_rate2[1:4,], col=alpha("forestgreen",0.5),border="Forestgreen",
        at=c(1:8)-0.3, 
        ylim=c(0,1),xlim=c(0,9),xaxt="n",cex=0.25,boxwex=0.25,xlab="Age (day)",ylab="Survial rate")
boxplot(suv_rate2[5:8,], add=T,col=alpha("salmon",0.5),border="salmon",
        at=1:8,axes=F,boxwex=0.25,cex=0.25)
boxplot(suv_rate2[9:12,], add=T,col=alpha("royalblue",0.5),border="royalblue",
        at=c(1:8)+0.3,axes=F,boxwex=0.25,cex=0.25)
axis(1,at=c(1:8),labels = c(1:8)*7)
#dev.off()

suv_rate3=read.csv("../surv_rate_new_cold.csv",header = T,row.names = 1,stringsAsFactors = F)

#png("../survival_rate_newCGE_cold.png",width = 8.7,height = 6,units = "cm",pointsize = 6,res=600)
par(mar=c(5,5,2,2))
boxplot(suv_rate3[1:4,], col=alpha("forestgreen",0.5),border="Forestgreen",
        at=c(1:8)-0.3, 
        ylim=c(0,1),xlim=c(0,9),xaxt="n",cex=0.25,boxwex=0.25,xlab="Age (day)",ylab="Survial rate")
boxplot(suv_rate3[5:8,], add=T,col=alpha("salmon",0.5),border="salmon",
        at=1:8,axes=F,boxwex=0.25,cex=0.25)
boxplot(suv_rate3[9:12,], add=T,col=alpha("royalblue",0.5),border="royalblue",
        at=c(1:8)+0.3,axes=F,boxwex=0.25,cex=0.25)
axis(1,at=c(1:8),labels = c(1:8)*12)
#dev.off()

####phenotype: life time fecundity####
fec_dat = read.csv("../final_eggs_all.csv")
fec_dat$day = as.numeric(as.factor(fec_dat$date))
fecPerDay = tapply(fec_dat$eggs,paste(fec_dat$day,fec_dat$pop,sep = "_"),sum,na.rm = T)
fecDay = as.numeric(strsplit2(names(fecPerDay),"_")[,1])
fecEvo = substr(strsplit2(names(fecPerDay),"_")[,2],1,3)
fecEvo[fecEvo=="rap"] = "anc"
fecSum = tapply(fec_dat$eggs,fec_dat$pop,sum,na.rm = T)
png("supplementary_life-time_fecundity_full.png",width = 17.4,height = 8.7,units = "cm",res = 600,pointsize = 8)
par(mar = c(4,4,1,1))
boxplot(fecPerDay~fecEvo+fecDay,border = c("gold","purple"),cex = 0.1,col=NA,xaxt = "n",xlab = "Age (Day)",ylab = "# of egg per population")
axis(1,at = seq(1.5,130.5,2),labels = sort(unique(fecDay))+1)
dev.off()

png("fig2_v2.png",width = 17.4,height = 17.4,units = "cm",res = 600,pointsize = 8)
layout(matrix(c(1,1,2,3),2,2,byrow = T))
layout.show(3)
par(mar=c(4,4,2,2))
boxplot(suv_rate[1:6,-1:-2], col=alpha("gold",0.5),border="gold",
        at=c(1:8,seq(10,18,2),seq(23,58,5))-0.25, 
        ylim=c(0,1),xlim=c(0,60),xaxt="n",cex=0.4,boxwex=0.4,xlab="Age (day)",ylab="Survial rate")
boxplot(suv_rate[-1:-6,-1:-2], add=T,col=alpha("purple",0.5),border="purple",
        at=c(1:8,seq(10,18,2),seq(23,58,5))+0.25,axes=F,boxwex=0.4,cex=0.4)
axis(1,at=c(1:8,seq(10,18,2),seq(23,58,5)))
par(mar=c(4,4,2,2))
boxplot(fecPerDay[fecDay<15]~fecEvo[fecDay<15]+fecDay[fecDay<15],
        border = c("gold","purple"),col=alpha(c("gold","purple"),0.5),
        cex = 0.4,xaxt = "n",xlab = "Age (Day)",ylab = "# of egg per population",lex.order =F)
axis(1,at = seq(1.5,28.5,2),labels = 2:15)
par(mar=c(4,4,2,2))
boxplot(fecSum~rep(c("Evolved","Ancestral"),each =5),
        border = c("gold","purple"),col=alpha(c("gold","purple"),0.5),
        xlab = "Population",ylab = "Life time fecunditiy (# of eggs)")
dev.off()

####gene expression data input and processing####
count_dat=read.csv("age_formal_fullset_readcounts.csv",header = T,row.names = 1)
count_dat_use=count_dat[,setdiff(colnames(count_dat),c("B02_12","B02_30","B03_1","H12_30"))]
colnames(count_dat_use)[grep("H12_30_t",colnames(count_dat_use))]="H12_30"
count_dat_use=count_dat_use[,-grep("_t",colnames(count_dat_use))]
count_dat_use=count_dat_use[,-grep("_60",colnames(count_dat_use))]
count_dat_filtered=count_dat_use[apply(cpm(count_dat_use),1,function(x) !sum(x<0.1)>=1),]
group=paste(substr(colnames(count_dat_filtered),1,1),substr(colnames(count_dat_filtered),5,6),sep="_")
y=DGEList(counts = count_dat_filtered,group = group)
y=calcNormFactors(y)
plotMDS(y)
#pca
pca=prcomp(t(log(cpm(y))))
ve=(pca$sdev^2)/sum(pca$sdev^2)
mycol=rainbow(60)[as.numeric(substr(group,3,4))+3]
mycol2=rainbow(60)[round(mort_rate_sample[colnames(count_dat_filtered),6]*50+3,digits = 0)]
#mycol2[is.na(mycol2)]="grey"
mypch=ifelse(substr(group,1,1)%in%"B",1,19)
#png("PCA.png",height = 8,width = 12,units = "cm",pointsize = 8,res = 600)
par(mar=c(5,5,2,2))
plot(pca$x,col=mycol,pch=mypch,asp=1,ylim=c(-10,35),
     xlab=paste0("PC1 (",round(ve[1],3)*100,"%)"),ylab=paste0("PC2 (",round(ve[2],3)*100,"%)"))
legend("bottomleft",legend = c("base","evolved"),pch=c(1,19),bty = "n")
legend("top",legend = sort(as.numeric(unique(substr(group,3,4)))),bty="n",
       fill=rainbow(60)[sort(as.numeric(unique(substr(group,3,4))))],horiz = T,x.intersp = 0.1)
#dev.off()
#png("PCA_1vs3.png",height = 12,width = 12,units = "cm",pointsize = 8,res = 600)
plot(pca$x[,c(1,3)],col=mycol,pch=mypch,asp=1,
     xlab=paste0("PC1 (",round(ve[1],3)*100,"%)"),ylab=paste0("PC3 (",round(ve[3],3)*100,"%)"))
legend("topleft",legend = c("base","evolved"),pch=c(1,19),bty = "n")
legend("bottom",legend = sort(as.numeric(unique(substr(group,3,4)))),bty="n",
       fill=rainbow(60)[sort(as.numeric(unique(substr(group,3,4))))],horiz = T,x.intersp = 0.1)
#dev.off()
#png("PCA_comparison.png",height = 12,width = 24,units = "cm",pointsize = 8,res = 600)
par(mfrow=c(1,2))
plot(pca$x,col=mycol,pch=mypch,asp=1,
     xlab=paste0("PC1 (",round(ve[1],3)*100,"%)"),ylab=paste0("PC2 (",round(ve[2],3)*100,"%)"))
legend("topleft",legend = c("base","evolved"),pch=c(1,19),bty = "n")
legend("bottom",legend = sort(as.numeric(unique(substr(group,3,4)))),bty="n",
       fill=rainbow(60)[sort(as.numeric(unique(substr(group,3,4))))+3],horiz = T,x.intersp = 0.1)
plot(pca$x,col=mycol2,pch=mypch,asp=1,
     xlab=paste0("PC1 (",round(ve[1],3)*100,"%)"),ylab=paste0("PC2 (",round(ve[2],3)*100,"%)"))
legend("topleft",legend = c("base","evolved"),pch=c(1,19),bty = "n")
legend("bottom",legend = seq(0,0.8,0.1),bty="n",
       fill=rainbow(60)[log(seq(0,0.8,0.1)+1)/log(1.8)*50+3],horiz = T,x.intersp = 0.1)
#dev.off()

png("PCA_comparison_1vs3.png",height = 12,width = 24,units = "cm",pointsize = 8,res = 600)
par(mfrow=c(1,2))
plot(pca$x[,c(1,3)],col=mycol,pch=mypch,asp=1,
     xlab=paste0("PC1 (",round(ve[1],3)*100,"%)"),ylab=paste0("PC3 (",round(ve[3],3)*100,"%)"))
legend("topleft",legend = c("base","evolved"),pch=c(1,19),bty = "n")
legend("bottom",legend = sort(as.numeric(unique(substr(group,3,4)))),bty="n",
       fill=rainbow(60)[sort(as.numeric(unique(substr(group,3,4))))+3],horiz = T,x.intersp = 0.1)
plot(pca$x[,c(1,3)],col=mycol2,pch=mypch,asp=1,
     xlab=paste0("PC1 (",round(ve[1],3)*100,"%)"),ylab=paste0("PC3 (",round(ve[3],3)*100,"%)"))
legend("topleft",legend = c("base","evolved"),pch=c(1,19),bty = "n")
legend("bottom",legend = seq(0,0.8,0.1),bty="n",
       fill=rainbow(60)[log(seq(0,0.8,0.1)+1)/log(1.8)*50+3],horiz = T,x.intersp = 0.1)
dev.off()

####edgeR analysis####
ModelDesign=model.matrix(~group)
DGE=estimateDisp(y,design = ModelDesign,robust = T)
GLM=glmFit(DGE,design = ModelDesign)
LRT_res=glmLRT(GLM,coef = 2:20)
res_table=LRT_res$table
res_table$padj=p.adjust(res_table$PValue,method = "BH")
sum(res_table$padj<0.05)
sig.gene=rownames(y)[res_table$padj<0.05]

ModelDesign2=model.matrix(~0+group)
DGE2=estimateDisp(y,design = ModelDesign2,robust = T)
GLM2=glmFit(DGE2,design = ModelDesign2)
mycontrast=makeContrasts(evo=(groupH_1+groupH_3+groupH_5+groupH_7+groupH_9+groupH_12+
                                  groupH_16+groupH_20+groupH_30+groupH_40+groupH_50)/11-
                           (groupB_1+groupB_3+groupB_5+groupB_7+groupB_9+groupB_12+
                              groupB_16+groupB_20+groupB_30+groupB_40+groupB_50)/11,
                         age1=(groupH_1+groupB_1)/2-(groupH_3+groupB_3)/2,
                         age2=(groupH_1+groupB_1)/2-(groupH_5+groupB_5)/2,
                         age3=(groupH_1+groupB_1)/2-(groupH_7+groupB_7)/2,
                         age4=(groupH_1+groupB_1)/2-(groupH_9+groupB_9)/2,
                         age5=(groupH_1+groupB_1)/2-(groupH_12+groupB_12)/2,
                         age6=(groupH_1+groupB_1)/2-(groupH_16+groupB_16)/2,
                         age7=(groupH_1+groupB_1)/2-(groupH_20+groupB_20)/2,
                         age8=(groupH_1+groupB_1)/2-(groupH_30+groupB_30)/2,
                         age9=(groupH_1+groupB_1)/2-(groupH_40+groupB_40)/2,
                         age10=(groupH_1+groupB_1)/2-(groupH_50+groupB_50)/2,
                         levels = ModelDesign2)
LRT_res_evo=glmLRT(GLM2,contrast = mycontrast[,"evo"])
res_table_evo=LRT_res_evo$table
res_table_evo$padj=p.adjust(res_table_evo$PValue,method = "BH")
sum(res_table_evo$padj<0.05)
LRT_res_age=glmLRT(GLM2,contrast = mycontrast[,2:11])
res_table_age=LRT_res_age$table
res_table_age$padj=p.adjust(res_table_age$PValue,method = "BH")
sum(res_table_age$padj<0.05)
venn.diagram(list(rownames(y)[res_table_age$padj<0.05],rownames(y)[res_table_evo$padj<0.05]),
             filename = "venn_linear_modeling.png",height = 12,width = 12,units = "cm",resolution = 600,
             imagetype = "png",category.names = c("Age: 10513","Evolution: 5439"),pointsize=8,cat.pos=0)

evo=strsplit2(group,"_")[,1]
age=strsplit2(group,"_")[,2]
ModelDesign3=model.matrix(~evo+age+evo:age)
DGE3=estimateDisp(y,design = ModelDesign3,robust = T)
GLM3=glmFit(DGE3,design = ModelDesign3)
LRT_res_evo2=glmLRT(GLM3,coef = 2)
res_table_evo2=LRT_res_evo2$table
res_table_evo2$padj=p.adjust(res_table_evo2$PValue,method = "BH")
sum(res_table_evo2$padj<0.05)

LRT_res_age2=glmLRT(GLM3,coef = 3:10)
res_table_age2=LRT_res_age2$table
res_table_age2$padj=p.adjust(res_table_age2$PValue,method = "BH")
sum(res_table_age2$padj<0.05)

LRT_res_inter=glmLRT(GLM3,coef = 11:20)
res_table_inter=LRT_res_inter$table
res_table_inter$padj=p.adjust(res_table_inter$PValue,method = "BH")
sum(res_table_inter$padj<0.05)

#not working as the design matrix is not full rank...
pseudo_age=as.factor(mort_rate_sample[colnames(y),]$sudo_age)
#levels(pseudo_age)=1:13
# group_pseudo_age=paste(evo,pseudo_age,sep = "_")
# y2=DGEList(counts = count_dat_filtered,group = group_pseudo_age)
# y2=calcNormFactors(y2)
# 
# ModelDesign4=model.matrix(~evo+pseudo_age+evo:pseudo_age)
# DGE4=estimateDisp(y2,design = ModelDesign4,robust = T)
# GLM4=glmFit(DGE4,design = ModelDesign4)
# LRT_res_evo3=glmLRT(GLM4,coef = 2)
# res_table_evo3=LRT_res_evo3$table
# res_table_evo3$padj=p.adjust(res_table_evo3$PValue,method = "BH")
# sum(res_table_evo3$padj<0.05)
# 
# LRT_res_age3=glmLRT(GLM4,coef = 3:10)
# res_table_age3=LRT_res_age3$table
# res_table_age3$padj=p.adjust(res_table_age3$PValue,method = "BH")
# sum(res_table_age3$padj<0.05)
# 
# LRT_res_inter2=glmLRT(GLM4,coef = 11:20)
# res_table_inter2=LRT_res_inter2$table
# res_table_inter2$padj=p.adjust(res_table_inter2$PValue,method = "BH")
# sum(res_table_inter2$padj<0.05)


####GO analysis on DE genes####
query_ID_evo=list(evo_up=rownames(y)[res_table_evo$padj<0.05&res_table_evo$logFC>0],
                  evo_dn=rownames(y)[res_table_evo$padj<0.05&res_table_evo$logFC<0])
GO_res_table=list()
for (i in names(query_ID_evo)){
  idx=rownames(y)%in%query_ID_evo[[i]]
  tmp=factor(as.integer(idx))
  names(tmp)=rownames(y)#genelist#
  tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")#data preparation#
  resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")#enrichment test#
  resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)#analysis of results#
  GO_res_table[[i]]=tmp_res
}

head(GO_res_table$evo_up)
head(GO_res_table$evo_dn,30)


####WGCNA_module_construction####
allowWGCNAThreads(24)
expr_dat=t(log(cpm(y)))
expr_dat_B=t(log(cpm(y)[,substr(group,1,1)%in%"B"]))
expr_dat_H=t(log(cpm(y)[,substr(group,1,1)%in%"H"]))

#B
powers = 1:20
sft_dat_B = pickSoftThreshold(expr_dat_B, powerVector = powers,verbose = 5)

plot(sft_dat_B$fitIndices[,1], -sign(sft_dat_B$fitIndices[,3])*sft_dat_B$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_dat_B$fitIndices[,1], -sign(sft_dat_B$fitIndices[,3])*sft_dat_B$fitIndices[,2],
     labels=powers,cex=1,col="red")
abline(h=0.90,col="red")

plot(sft_dat_B$fitIndices[,1], sft_dat_B$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_dat_B$fitIndices[,1], sft_dat_B$fitIndices[,5], labels=powers, cex=1,col="red")

net_dat_B = blockwiseModules(expr_dat_B, power = 4,maxBlockSize = ncol(expr_dat_B),
                           TOMType = "signed",networkType = "signed", minModuleSize = 20,
                           reassignThreshold = 1e-6, mergeCutHeight = 0.15,
                           numericLabels = TRUE, pamRespectsDendro = T,
                           saveTOMs = T,nThreads = 12,
                           verbose = 3)

table(net_dat_B$colors)
mergedColors_B = labels2colors(net_dat_B$colors)

plotDendroAndColors(net_dat_B$dendrograms[[1]], mergedColors_B[net_dat_B$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

load("./blockwiseTOM-block.1.RData")
plotTOM = as.matrix((1-TOM)^4)
diag(plotTOM) = NA
png("./heatmap_TOM_WGCNA_B.png",height = 16,width = 16,units = "cm",res = 600)
TOMplot(plotTOM, net_dat_B$dendrograms[[1]], mergedColors_B[net_dat_B$blockGenes[[1]]],
        main = "Network heatmap plot, all genes",col=heat.colors(100))
dev.off()

set.seed(100)
eg.idx=sample(which(net_dat_B$colors==2),2)
png("example_coexpression.png",width = 8.7,height = 8.7,units = "cm",pointsize = 8,res=600)
plot(expr_dat_B[,eg.idx],asp=1,col=mycol[substr(group,1,1)%in%"B"])
dev.off()
cor(expr_dat_B[,eg.idx])

#H
powers = 1:20
sft_dat_H = pickSoftThreshold(expr_dat_H, powerVector = powers, verbose = 5)

plot(sft_dat_H$fitIndices[,1], -sign(sft_dat_H$fitIndices[,3])*sft_dat_H$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_dat_H$fitIndices[,1], -sign(sft_dat_H$fitIndices[,3])*sft_dat_H$fitIndices[,2],
     labels=powers,cex=1,col="red")
abline(h=0.90,col="red")

plot(sft_dat_H$fitIndices[,1], sft_dat_H$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_dat_H$fitIndices[,1], sft_dat_H$fitIndices[,5], labels=powers, cex=1,col="red")

net_dat_H = blockwiseModules(expr_dat_H, power = 4,maxBlockSize = ncol(expr_dat_H),
                             TOMType = "signed",networkType = "signed", minModuleSize = 20,
                             reassignThreshold = 1e-6, mergeCutHeight = 0.15,
                             numericLabels = TRUE, pamRespectsDendro = T,
                             saveTOMs = F,nThreads = 12,
                             verbose = 3)

table(net_dat_H$colors)
mergedColors_H = labels2colors(net_dat_H$colors)

plotDendroAndColors(net_dat_H$dendrograms[[1]], mergedColors_H[net_dat_H$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#overall
powers = 1:20
sft_dat = pickSoftThreshold(expr_dat, powerVector = powers, verbose = 5)

plot(sft_dat$fitIndices[,1], -sign(sft_dat$fitIndices[,3])*sft_dat$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_dat$fitIndices[,1], -sign(sft_dat$fitIndices[,3])*sft_dat$fitIndices[,2],
     labels=powers,cex=1,col="red")
abline(h=0.90,col="red")

plot(sft_dat$fitIndices[,1], sft_dat$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_dat$fitIndices[,1], sft_dat$fitIndices[,5], labels=powers, cex=1,col="red")

net_dat = blockwiseModules(expr_dat, power = 3,maxBlockSize = ncol(expr_dat),
                             TOMType = "signed",networkType = "signed", minModuleSize = 20,
                             reassignThreshold = 1e-5, mergeCutHeight = 0.10,
                             numericLabels = TRUE, pamRespectsDendro = T,
                             saveTOMs = F,nThreads = 12,
                             verbose = 3)

table(net_dat$colors)
mergedColors = labels2colors(net_dat$colors)

plotDendroAndColors(net_dat$dendrograms[[1]], mergedColors[net_dat$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



age_B=substr(group,3,4)[substr(group,1,1)%in%"B"]
age_H=substr(group,3,4)[substr(group,1,1)%in%"H"]

scl_expr_dat=apply(expr_dat,2,scale)
scl_avg_expr_dat=apply(scl_expr_dat,2,function(x) tapply(x,group,mean))

scl_expr_dat_B=apply(expr_dat_B,2,scale)
scl_avg_expr_dat_B=apply(scl_expr_dat_B,2,function(x) tapply(x,age_B,mean))
scl_avg_expr_dat_B=scl_avg_expr_dat_B[as.character(sort(as.numeric(rownames(scl_avg_expr_dat_B)))),]
scl_expr_dat_H=apply(expr_dat_H,2,scale)
scl_avg_expr_dat_H=apply(scl_expr_dat_H,2,function(x) tapply(x,age_H,mean))
scl_avg_expr_dat_H=scl_avg_expr_dat_H[as.character(sort(as.numeric(rownames(scl_avg_expr_dat_H)))),]

scl_avg_expr_dat_BB=scl_avg_expr_dat[grep("B",rownames(scl_avg_expr_dat)),][paste0("B_",sort(as.numeric(rownames(scl_avg_expr_dat_B)))),]
scl_avg_expr_dat_HH=scl_avg_expr_dat[grep("H",rownames(scl_avg_expr_dat)),][paste0("H_",sort(as.numeric(rownames(scl_avg_expr_dat_H)))),]

png("./WGCNA_module_expression_B.png",width = 17.4,height = 8.7,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,4))
for (i in 0:7){
  idx=net_dat_B$colors==i
  plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Module",i))
  axis(1,at=sort(as.numeric(unique(age_B))))
  apply(scl_avg_expr_dat_B[,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("grey20",0.2)))
}
dev.off()

png("./WGCNA_module_expression_B_range.png",width = 17.4,height = 8.7,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,4))
for (i in 0:7){
  idx=net_dat_B$colors==i
  avg_B=apply(scl_avg_expr_dat_B[,idx],1,mean)
  sd_B=apply(scl_avg_expr_dat_B[,idx],1,sd)

  plot(NA,xlim=c(0,50),ylim=c(-1.5,1.5),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Module",i))
  axis(1,at=sort(as.numeric(unique(age_H))))
  lines(sort(as.numeric(unique(age_B))),avg_B,col="forestgreen",lwd=3)
  polygon(x=c(sort(as.numeric(unique(age_B))),rev(sort(as.numeric(unique(age_B))))),
          y = c(avg_B-sd_B,rev(avg_B+sd_B)),col=alpha("forestgreen",0.3),border = "forestgreen",lty=2)
}
dev.off()

png("./WGCNA_module_expression_H.png",width = 13.05,height = 8.7,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,3))
for (i in 0:5){
  idx=net_dat_H$colors==i
  plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Module",i))
  axis(1,at=sort(as.numeric(unique(age_H))))
  apply(scl_avg_expr_dat_H[,idx],2,function(x) lines(sort(as.numeric(unique(age_H))),x,col=alpha("grey20",0.2)))
}
dev.off()

png("./WGCNA_module_expression_H_range.png",width = 13.05,height = 8.7,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,3))
for (i in 0:5){
  idx=net_dat_H$colors==i
  avg_H=apply(scl_avg_expr_dat_H[,idx],1,mean)
  sd_H=apply(scl_avg_expr_dat_H[,idx],1,sd)
  
  plot(NA,xlim=c(0,50),ylim=c(-1.5,1.5),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Module",i))
  axis(1,at=sort(as.numeric(unique(age_H))))
  lines(sort(as.numeric(unique(age_B))),avg_H,col="salmon",lwd=3)
  polygon(x=c(sort(as.numeric(unique(age_B))),rev(sort(as.numeric(unique(age_B))))),
          y = c(avg_H-sd_H,rev(avg_H+sd_H)),col=alpha("salmon",0.3),border = "salmon",lty=2)
}
dev.off()

png("./WGCNA_module_expression_combined.png",width = 17.4,height = 8.7,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,4))
for (i in 0:7){
  idx=net_dat$colors==i
  avg_BB=apply(scl_avg_expr_dat_BB[,idx],1,mean)
  sd_BB=apply(scl_avg_expr_dat_BB[,idx],1,sd)
  avg_HH=apply(scl_avg_expr_dat_HH[,idx],1,mean)
  sd_HH=apply(scl_avg_expr_dat_HH[,idx],1,sd)
  
  plot(NA,xlim=c(0,50),ylim=c(-1.5,1.5),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Module",i))
  axis(1,at=sort(as.numeric(unique(age_H))))
  lines(sort(as.numeric(unique(age_B))),avg_BB,col="forestgreen",lwd=3)
  polygon(x=c(sort(as.numeric(unique(age_B))),rev(sort(as.numeric(unique(age_B))))),
          y = c(avg_BB-sd_BB,rev(avg_BB+sd_BB)),col=alpha("forestgreen",0.3),border = "forestgreen",lty=2)
  lines(sort(as.numeric(unique(age_H))),avg_HH,col="salmon",lwd=3)
  polygon(x=c(sort(as.numeric(unique(age_H))),rev(sort(as.numeric(unique(age_H))))),
          y = c(avg_HH-sd_HH,rev(avg_HH+sd_HH)),col=alpha("salmon",0.3),border = "salmon",lty = 2)
}
dev.off()

png("./WGCNA_module_expression_comparison_B_net_based.png",width = 17.4,height = 8.7,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,4))
for (i in 0:7){
  idx=net_dat_B$colors==i
  avg_BB=apply(scl_avg_expr_dat_BB[,idx],1,mean)
  sd_BB=apply(scl_avg_expr_dat_BB[,idx],1,sd)
  avg_HH=apply(scl_avg_expr_dat_HH[,idx],1,mean)
  sd_HH=apply(scl_avg_expr_dat_HH[,idx],1,sd)
  
  plot(NA,xlim=c(0,50),ylim=c(-1.5,1.5),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Module",i))
  axis(1,at=sort(as.numeric(unique(age_H))))
  lines(sort(as.numeric(unique(age_B))),avg_BB,col="forestgreen",lwd=3)
  polygon(x=c(sort(as.numeric(unique(age_B))),rev(sort(as.numeric(unique(age_B))))),
          y = c(avg_BB-sd_BB,rev(avg_BB+sd_BB)),col=alpha("forestgreen",0.3),border = "forestgreen",lty=2)
  lines(sort(as.numeric(unique(age_H))),avg_HH,col="salmon",lwd=3)
  polygon(x=c(sort(as.numeric(unique(age_H))),rev(sort(as.numeric(unique(age_H))))),
          y = c(avg_HH-sd_HH,rev(avg_HH+sd_HH)),col=alpha("salmon",0.3),border = "salmon",lty = 2)
}
dev.off()

png("./WGCNA_module_expression_comparison_B_net_based_rescaled.png",width = 17.4,height = 8.7,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,4))
for (i in 0:7){
  idx=net_dat_B$colors==i
  avg_BB=apply(scl_avg_expr_dat_BB[,idx],1,mean)
  sd_BB=apply(scl_avg_expr_dat_BB[,idx],1,sd)
  avg_HH=apply(scl_avg_expr_dat_HH[,idx],1,mean)
  sd_HH=apply(scl_avg_expr_dat_HH[,idx],1,sd)
  
  plot(NA,xlim=c(0,0.6),ylim=c(-1.5,1.5),xlab="Mortality rate",ylab="normalized mean expression",xaxt="n",
       main=paste0("Module",i))
  axis(1)
  lines(avg_mort_rate_sample[1,-12],avg_BB,col="forestgreen",lwd=3)
  polygon(x=c(avg_mort_rate_sample[1,-12],rev(avg_mort_rate_sample[1,-12])),
          y = c(avg_BB-sd_BB,rev(avg_BB+sd_BB)),col=alpha("forestgreen",0.3),border = "forestgreen",lty=2)
  lines(avg_mort_rate_sample[2,-12],avg_HH,col="salmon",lwd=3)
  polygon(x=c(avg_mort_rate_sample[2,-12],rev(avg_mort_rate_sample[2,-12])),
          y = c(avg_HH-sd_HH,rev(avg_HH+sd_HH)),col=alpha("salmon",0.3),border = "salmon",lty = 2)
}
dev.off()


####WGCNA_module GO analysis####
mod_dat_GO_res=list()
for (i in sort(unique(net_dat$colors))){
  idx=net_dat$colors==i
  tmp=factor(as.integer(idx))
  names(tmp)=rownames(y)#genelist#
  tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")#data preparation#
  resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")#enrichment test#
  resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)#analysis of results#
  mod_dat_GO_res[[i+1]]=tmp_res
}

names(mod_dat_GO_res)=paste0("Module",0:7)
lapply(mod_dat_GO_res,function(x) x[1:3,c(1:2,8)])

for (i in 1:length(mod_dat_GO_res)){
  write.table(mod_dat_GO_res[[i]],paste0("./WGCNA_GO/",names(mod_dat_GO_res)[i],"_GO_res.txt"),quote = F,
              row.names = F,sep="\t")
}
####WGCNA_module_comparison####
multiexpr=list(list(data=expr_dat_B),list(data=expr_dat_H))
modulelist=list(mergedColors_B,mergedColors_B)
MP=modulePreservation(multiexpr, modulelist,
                      referenceNetworks = c(1:2),
                      loadPermutedStatistics = FALSE,
                      nPermutations = 100,
                      networkType = "signed",
                      maxModuleSize = 10000,
                      maxGoldModuleSize = 10000,
                      verbose = 3)

# Calculate the contingency table and p-values
overlap = overlapTable(mergedColors_B, mergedColors_H)
# The numMat will encode color. We use -log of the p value.
numMat = -log10(overlap$pTable)
numMat[numMat >50] = 50
# Prepare for generating a color-coded plot of the overlap table. The text of the table will consist of
# counts and corresponding p-values.
textMat = paste(overlap$countTable, "\n", signif(overlap$pTable, 2));
dim(textMat) = dim(numMat)
# Additional information for the plot. These will be used shortly.
xLabels = paste("M", sort(unique(mergedColors_H)));
yLabels = paste("M", sort(unique(mergedColors_B)));
xSymbols = paste(sort(unique(mergedColors_H)), ": ", table(mergedColors_H), sep = "")
ySymbols = paste(sort(unique(mergedColors_B)), ": ", table(mergedColors_B), sep = "")

png("module_preservation.png", w = 17.4, h = 17.4,units = "cm",res = 600,pointsize = 10)
fp = TRUE
layout(matrix(c(1,2,5, 3,4,5), 3, 2),
       heights = c(3, 1, 6), widths = c(1, 1));
#layout.show(5);
par(mgp = c(3, 1, 0));
plotDendroAndColors(net_dat_B$dendrograms[[1]],
                    cbind(mergedColors_B, mergedColors_H),
                    c("Ancestral modules", "Evolved modules"),
                    setLayout = FALSE,
                    marAll = c(1, 6, 2.7, 0.2),
                    addGuide = FALSE,
                    main = "Ancestral gene dendrogram\nand module preservation", cex.main = 1.2,
                    dendroLabels = FALSE, hang = 0.03, cex.colorLabels = 0.7, abHeight = 0.95)
par(mgp = c(3, 1, 0));
plotDendroAndColors(net_dat_H$dendrograms[[1]],
                    cbind(mergedColors_B, mergedColors_H),
                    c("Ancestral modules", "Evolved modules"),
                    setLayout = FALSE,
                    marAll = c(1, 6, 2.7, 0.2),
                    addGuide = FALSE,
                    main = "Evolved gene dendrogram\nand module preservation", cex.main = 1.2,
                    dendroLabels = FALSE, hang = 0.03, cex.colorLabels = 0.7, abHeight = 0.95)

# Plot the overlap table
fcex = 1.00;
pcex = 1.0
fcexl = 1.00;
pcexl = 1.00;
par(mar = c(8,8, 2, 1.0));
labeledHeatmap(Matrix = numMat,
               xLabels = xLabels, xSymbols = xSymbols,
               yLabels = yLabels, ySymbols = ySymbols,
               colorLabels = TRUE,
               colors = greenWhiteRed(100)[50:100],
               textMatrix = textMat, cex.text = if (fp) fcex else pcex, setStdMargins = FALSE,
               cex.lab = if (fp) fcexl else pcexl,
               xColorWidth = 0.08,
               main = "Ancestral modules (rows) vs. Evolved modules (columns)", cex.main = 1.2)
dev.off()

####maSigPro####
expr_design=data.frame(age=as.numeric(substr(colnames(y),5,6)),
                       replicate=as.numeric(as.factor(group)),
                       anc=as.numeric(substr(colnames(y),1,1)=="B"),
                       evo=as.numeric(substr(colnames(y),1,1)=="H"))
rownames(expr_design)=colnames(y)

masigpro.design = make.design.matrix(expr_design, degree = 2)
fit.logcpm = p.vector(t(expr_dat), masigpro.design, Q = 0.05, MT.adjust = "BH", min.obs = 20)
# fit.cpm = p.vector(cpm(y), masigpro.design, Q = 0.05, MT.adjust = "BH", min.obs = 20,counts = T)

tstep.logcpm = T.fit(fit.logcpm, step.method = "backward", alfa = 0.05)
# tstep.cpm = T.fit(fit.cpm, step.method = "backward", alfa = 0.05)

sigs.logcpm = get.siggenes(tstep.logcpm, rsq = 0.6, vars = "groups")
# sigs.cpm = get.siggenes(tstep.cpm, rsq = 0.6, vars = "groups")

#see.genes(sigs.logcpm$sig.genes$evovsanc, show.fit = T, dis =expr_design$dis,
#          cluster.method="hclust" ,cluster.data = 1, k = 9) #issue with X11()...

sigset=apply(sigs.logcpm$sig.genes$evovsanc$sig.pvalues[,-1:-3],2,function(x) return(x<0.05))
sigset[is.na(sigset)]=0
sigset=as.data.frame(sigset)
colnames(sigset)=c("evo","age","agexevo","age2","age2xevo")
png("maSigPro_significance.png",width = 17.4,height = 8.7,units = "cm",res = 600,pointsize = 8)
UpSetR::upset(sigset,sets = c("evo","age","agexevo","age2","age2xevo"),
      group.by = "degree",keep.order = T)
dev.off()

png("maSigPro_significance_subset.png",width = 17.4,height = 8.7,units = "cm",res = 600,pointsize = 8)
UpSetR::upset(sigset,sets = c("evo","agexevo","age2xevo"),
      group.by = "degree",keep.order = T)
dev.off()

interIdx=apply(sigset[,c(3,5)],1,function(x) any(as.logical(x)))

EvoIdx=c()
for (i in 1:nrow(sigset)) {
  if (sigset[i,1]==1&sum(sigset[i,])==1) {
    EvoIdx[i]=TRUE
  }
  else EvoIdx[i]=FALSE
}


clusterdata=sigs.logcpm$sig.genes$evovsanc$sig.profiles[interIdx,]
dcorrel = matrix(rep(1, nrow(clusterdata)^2), nrow(clusterdata), 
                 nrow(clusterdata)) - cor(t(clusterdata),use = "pairwise.complete.obs")
outTab = sigs.logcpm$sig.genes$evovsanc$sig.pvalues
write.table(cbind(outTab[interIdx,],cut),
            "~/Dropbox (PopGen)/Wei-Yun/project_aging/TableS1b.txt",
            sep = "\t",row.names = T,col.names = T,quote = F)

#clusterdata_noninter=sigs.logcpm$sig.genes$evovsanc$sig.profiles[!interIdx,] #wylai add a new dataset not including the interaction term [2023.06.28]
#dcorrel_noninter = matrix(rep(1, nrow(clusterdata_noninter)^2), nrow(clusterdata_noninter), 
#nrow(clusterdata_noninter)) - cor(t(clusterdata_noninter),use = "pairwise.complete.obs")

#pheatmap(dcorrel) #2765*2765
clust = hclust(as.dist(dcorrel), method = "ward.D2")
cut = cutree(clust, k = 12)
cut[c("FBgn0035995","FBgn0040022","FBgn0250848","FBgn0030509","FBgn0030509")]

#clust_noninter = hclust(as.dist(dcorrel_noninter), method = "ward.D2")
#cut_noninter = cutree(clust_noninter, k = 4)

# png("./maSigPro_v2.png",width = 17.4,height = 17.4,units = "cm",res = 600,pointsize = 8)
# par(mfrow=c(3,4))
# for (i in sort(unique(cut))){
#   idx=cut==c(1,2,6,11,3,7,4,5,10,8,9,12)[i]
#   plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
#        main=paste0("Cluster",i,": ",sum(idx)))
#   axis(1,at=sort(as.numeric(unique(age_B))))
#   apply(scl_avg_expr_dat_BB[,rownames(clusterdata)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("green",0.2)))
#   apply(scl_avg_expr_dat_HH[,rownames(clusterdata)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("salmon",0.2)))
# }
# dev.off()

png("./maSigPro_ranged_v3.png",width = 17.4,height = 8.7*1.5,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(3,4))
for (i in sort(unique(cut))){
  idx=cut== c(1,2,6,11,3,7,4,5,10,8,9,12)[i]
  avg_BB=apply(scl_avg_expr_dat_BB[,rownames(clusterdata)][,idx],1,mean)
  sd_BB=apply(scl_avg_expr_dat_BB[,rownames(clusterdata)][,idx],1,sd)
  avg_HH=apply(scl_avg_expr_dat_HH[,rownames(clusterdata)][,idx],1,mean)
  sd_HH=apply(scl_avg_expr_dat_HH[,rownames(clusterdata)][,idx],1,sd)
  
  plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Cluster",i,": ",sum(idx)))
  axis(1,at=sort(as.numeric(unique(age_B))))
  lines(sort(as.numeric(unique(age_B))),avg_BB,col="gold",lwd=3)
  polygon(x=c(sort(as.numeric(unique(age_B))),rev(sort(as.numeric(unique(age_B))))),
          y = c(avg_BB-sd_BB,rev(avg_BB+sd_BB)),col=alpha("gold",0.3),border = "gold",lty=2)
  lines(sort(as.numeric(unique(age_H))),avg_HH,col="purple",lwd=3)
  polygon(x=c(sort(as.numeric(unique(age_H))),rev(sort(as.numeric(unique(age_H))))),
          y = c(avg_HH-sd_HH,rev(avg_HH+sd_HH)),col=alpha("purple",0.3),border = "purple",lty = 2)
}
dev.off()


png("./maSigPro_ranged_noninter.png",width = 17.4,height = 8.7*1.5,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,2))
for (i in sort(unique(cut_noninter))){
  idx=cut_noninter== i
  avg_BB=apply(scl_avg_expr_dat_BB[,rownames(clusterdata_noninter)][,idx],1,mean)
  sd_BB=apply(scl_avg_expr_dat_BB[,rownames(clusterdata_noninter)][,idx],1,sd)
  avg_HH=apply(scl_avg_expr_dat_HH[,rownames(clusterdata_noninter)][,idx],1,mean)
  sd_HH=apply(scl_avg_expr_dat_HH[,rownames(clusterdata_noninter)][,idx],1,sd)
  
  plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Cluster",i,": ",sum(idx)))
  axis(1,at=sort(as.numeric(unique(age_B))))
  lines(sort(as.numeric(unique(age_B))),avg_BB,col="forestgreen",lwd=3)
  polygon(x=c(sort(as.numeric(unique(age_B))),rev(sort(as.numeric(unique(age_B))))),
          y = c(avg_BB-sd_BB,rev(avg_BB+sd_BB)),col=alpha("forestgreen",0.3),border = "forestgreen",lty=2)
  lines(sort(as.numeric(unique(age_H))),avg_HH,col="salmon",lwd=3)
  polygon(x=c(sort(as.numeric(unique(age_H))),rev(sort(as.numeric(unique(age_H))))),
          y = c(avg_HH-sd_HH,rev(avg_HH+sd_HH)),col=alpha("salmon",0.3),border = "salmon",lty = 2)
}
dev.off()





clusterdata1=sigs.logcpm$sig.genes$evovsanc[[1]][sigset$evo==1&sigset$agexevo==0&sigset$age2xevo==0,]
dcorrel1 = matrix(rep(1, nrow(clusterdata1)^2), nrow(clusterdata1), 
                  nrow(clusterdata1)) - cor(t(clusterdata1),use = "pairwise.complete.obs")
clust1 = hclust(as.dist(dcorrel1), method = "ward.D")
cut1 = cutree(clust1, k = 6)

png("./maSigPro_clust1.png",width = 17.4,height = 17.4,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,3))
for (i in sort(unique(cut1))){
  idx=cut1==i
  plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Cluster",i,": ",sum(idx)))
  axis(1,at=sort(as.numeric(unique(age_B))))
  apply(scl_avg_expr_dat_BB[,rownames(clusterdata1)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("green",0.2)))
  apply(scl_avg_expr_dat_HH[,rownames(clusterdata1)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("salmon",0.2)))
}
dev.off()

clusterdata2=sigs.logcpm$sig.genes$evovsanc[[1]][sigset$evo==0&sigset$agexevo==1&sigset$age2xevo==0,]
dcorrel2 = matrix(rep(1, nrow(clusterdata2)^2), nrow(clusterdata2), 
                  nrow(clusterdata2)) - cor(t(clusterdata2),use = "pairwise.complete.obs")
clust2 = hclust(as.dist(dcorrel2), method = "ward.D")
cut2 = cutree(clust2, k = 6)

png("./maSigPro_clust2.png",width = 17.4,height = 17.4,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,3))
for (i in sort(unique(cut2))){
  idx=cut2==i
  plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Cluster",i,": ",sum(idx)))
  axis(1,at=sort(as.numeric(unique(age_B))))
  apply(scl_avg_expr_dat_BB[,rownames(clusterdata2)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("green",0.2)))
  apply(scl_avg_expr_dat_HH[,rownames(clusterdata2)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("salmon",0.2)))
}
dev.off()

clusterdata3=sigs.logcpm$sig.genes$evovsanc[[1]][sigset$evo==0&sigset$agexevo==0&sigset$age2xevo==1,]
dcorrel3 = matrix(rep(1, nrow(clusterdata3)^2), nrow(clusterdata3), 
                  nrow(clusterdata3)) - cor(t(clusterdata3),use = "pairwise.complete.obs")
clust3 = hclust(as.dist(dcorrel3), method = "ward.D")
cut3 = cutree(clust3, k = 6)

png("./maSigPro_clust3.png",width = 17.4,height = 17.4,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,3))
for (i in sort(unique(cut3))){
  idx=cut3==i
  plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Cluster",i,": ",sum(idx)))
  axis(1,at=sort(as.numeric(unique(age_B))))
  apply(scl_avg_expr_dat_BB[,rownames(clusterdata3)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("green",0.2)))
  apply(scl_avg_expr_dat_HH[,rownames(clusterdata3)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("salmon",0.2)))
}
dev.off()

clusterdata4=sigs.logcpm$sig.genes$evovsanc[[1]][sigset$evo==1&sigset$agexevo==1&sigset$age2xevo==0,]
dcorrel4 = matrix(rep(1, nrow(clusterdata4)^2), nrow(clusterdata4), 
                  nrow(clusterdata4)) - cor(t(clusterdata4),use = "pairwise.complete.obs")
clust4 = hclust(as.dist(dcorrel4), method = "ward.D")
cut4 = cutree(clust4, k = 8)

png("./maSigPro_clust4.png",width = 17.4,height = 17.4,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,4))
for (i in sort(unique(cut4))){
  idx=cut4==i
  plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Cluster",i,": ",sum(idx)))
  axis(1,at=sort(as.numeric(unique(age_B))))
  apply(scl_avg_expr_dat_BB[,rownames(clusterdata4)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("green",0.2)))
  apply(scl_avg_expr_dat_HH[,rownames(clusterdata4)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("salmon",0.2)))
}
dev.off()

clusterdata5=sigs.logcpm$sig.genes$evovsanc[[1]][sigset$evo==1&sigset$agexevo==0&sigset$age2xevo==1,]
dcorrel5 = matrix(rep(1, nrow(clusterdata5)^2), nrow(clusterdata5), 
                  nrow(clusterdata5)) - cor(t(clusterdata5),use = "pairwise.complete.obs")
clust5 = hclust(as.dist(dcorrel5), method = "ward.D")
cut5 = cutree(clust5, k = 8)

png("./maSigPro_clust5.png",width = 17.4,height = 17.4,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,4))
for (i in sort(unique(cut5))){
  idx=cut5==i
  plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Cluster",i,": ",sum(idx)))
  axis(1,at=sort(as.numeric(unique(age_B))))
  apply(scl_avg_expr_dat_BB[,rownames(clusterdata5)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("green",0.2)))
  apply(scl_avg_expr_dat_HH[,rownames(clusterdata5)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("salmon",0.2)))
}
dev.off()

clusterdata6=sigs.logcpm$sig.genes$evovsanc[[1]][sigset$evo==0&sigset$agexevo==1&sigset$age2xevo==1,]
dcorrel6 = matrix(rep(1, nrow(clusterdata6)^2), nrow(clusterdata6), 
                  nrow(clusterdata6)) - cor(t(clusterdata6),use = "pairwise.complete.obs")
clust6 = hclust(as.dist(dcorrel6), method = "ward.D")
cut6 = cutree(clust6, k = 6)

png("./maSigPro_clust6.png",width = 17.4,height = 17.4,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,3))
for (i in sort(unique(cut6))){
  idx=cut6==i
  plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Cluster",i,": ",sum(idx)))
  axis(1,at=sort(as.numeric(unique(age_B))))
  apply(scl_avg_expr_dat_BB[,rownames(clusterdata6)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("green",0.2)))
  apply(scl_avg_expr_dat_HH[,rownames(clusterdata6)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("salmon",0.2)))
}
dev.off()

clusterdata7=sigs.logcpm$sig.genes$evovsanc[[1]][sigset$evo==1&sigset$agexevo==1&sigset$age2xevo==1,]
dcorrel7 = matrix(rep(1, nrow(clusterdata7)^2), nrow(clusterdata7), 
                  nrow(clusterdata7)) - cor(t(clusterdata7),use = "pairwise.complete.obs")
clust7 = hclust(as.dist(dcorrel7), method = "ward.D")
cut7 = cutree(clust7, k = 6)

png("./maSigPro_clust7.png",width = 17.4,height = 17.4,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,3))
for (i in sort(unique(cut7))){
  idx=cut7==i
  plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Cluster",i,": ",sum(idx)))
  axis(1,at=sort(as.numeric(unique(age_B))))
  apply(scl_avg_expr_dat_BB[,rownames(clusterdata7)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("green",0.2)))
  apply(scl_avg_expr_dat_HH[,rownames(clusterdata7)][,idx],2,function(x) lines(sort(as.numeric(unique(age_B))),x,col=alpha("salmon",0.2)))
}
dev.off()






####maSigPro_pseudo_age####
expr_design=data.frame(pseudo_age=as.numeric(pseudo_age),
                       replicate=as.numeric(as.factor(group)),
                       anc=as.numeric(substr(colnames(y),1,1)=="B"),
                       evo=as.numeric(substr(colnames(y),1,1)=="H"))
rownames(expr_design)=colnames(y)

masigpro.design = make.design.matrix(expr_design, degree = 2)
fit.logcpm = p.vector(t(expr_dat), masigpro.design, Q = 0.05, MT.adjust = "BH", min.obs = 20)
fit.cpm = p.vector(cpm(y), masigpro.design, Q = 0.05, MT.adjust = "BH", min.obs = 20,counts = T)

tstep.logcpm = T.fit(fit.logcpm, step.method = "backward", alfa = 0.05)
tstep.cpm = T.fit(fit.cpm, step.method = "backward", alfa = 0.05)

sigs.logcpm = get.siggenes(tstep.logcpm, rsq = 0.6, vars = "groups")
sigs.cpm = get.siggenes(tstep.cpm, rsq = 0.6, vars = "groups")

#see.genes(sigs.logcpm$sig.genes$evovsanc, show.fit = T, dis =expr_design$dis,
#          cluster.method="hclust" ,cluster.data = 1, k = 9) #issue with X11()...

clusterdata=sigs.logcpm$sig.genes$evovsanc[[1]]
dcorrel = matrix(rep(1, nrow(clusterdata)^2), nrow(clusterdata), 
                 nrow(clusterdata)) - cor(t(clusterdata),use = "pairwise.complete.obs")
clust = hclust(as.dist(dcorrel), method = "ward.D")
cut = cutree(clust, k = 4)

png("./maSigPro_rescale.png",width = 17.4,height = 17.4,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(2,2))
for (i in sort(unique(cut))){
  idx=cut==i
  plot(NA,xlim=c(0,11),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Cluster",i,": ",sum(idx)))
  axis(1)
  apply(scl_avg_expr_dat_BB[,rownames(clusterdata)][,idx],2,function(x) lines(round(tapply(as.numeric(pseudo_age),group,mean))[c(1,5,8,10,11,2,3,4,6,7,9)],x,col=alpha("green",0.2)))
  apply(scl_avg_expr_dat_HH[,rownames(clusterdata)][,idx],2,function(x) lines(round(tapply(as.numeric(pseudo_age),group,mean))[c(1,5,8,10,11,2,3,4,6,7,9)+11],x,col=alpha("salmon",0.2)))
}
dev.off()

png("./maSigPro_ranged.png",width = 17.4,height = 17.4,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(3,4))
for (i in sort(unique(cut))){
  idx=cut==i
  avg_BB=apply(scl_avg_expr_dat_BB[,rownames(clusterdata)][,idx],1,mean)
  sd_BB=apply(scl_avg_expr_dat_BB[,rownames(clusterdata)][,idx],1,sd)
  avg_HH=apply(scl_avg_expr_dat_HH[,rownames(clusterdata)][,idx],1,mean)
  sd_HH=apply(scl_avg_expr_dat_HH[,rownames(clusterdata)][,idx],1,sd)
  
  plot(NA,xlim=c(0,50),ylim=c(-3,3),xlab="age",ylab="normalized mean expression",xaxt="n",
       main=paste0("Cluster",i,": ",sum(idx)))
  axis(1,at=sort(as.numeric(unique(age_B))))
  lines(sort(as.numeric(unique(age_B))),avg_BB,col="forestgreen",lwd=3)
  polygon(x=c(sort(as.numeric(unique(age_B))),rev(sort(as.numeric(unique(age_B))))),
          y = c(avg_BB-sd_BB,rev(avg_BB+sd_BB)),col=alpha("forestgreen",0.3),border = "forestgreen",lty=2)
  lines(sort(as.numeric(unique(age_H))),avg_HH,col="salmon",lwd=3)
  polygon(x=c(sort(as.numeric(unique(age_H))),rev(sort(as.numeric(unique(age_H))))),
          y = c(avg_HH-sd_HH,rev(avg_HH+sd_HH)),col=alpha("salmon",0.3),border = "salmon",lty = 2)
}
dev.off()



####maSigPro_cluster GO analysis####
clu_ID=lapply(c(1,2,6,11,3,7,4,5,10,8,9,12),function(x) rownames(clusterdata)[cut==x])
clu_dat_GO_res=lapply(clu_ID,function(x) {
  tmp=factor(as.integer(rownames(y)%in%x))
  names(tmp)=rownames(y)#genelist#
  tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")#data preparation#
  resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")#enrichment test#
  resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)#analysis of results#
  return(tmp_res)
})

names(clu_dat_GO_res)=paste0("Cluster",1:12)
lapply(clu_dat_GO_res,function(x) x[1:10,c(1:2,8)])

dir.create("~/Dropbox (PopGen)/Wei-Yun/project_aging/maSigPro_GO_v2/")
for (i in 1:length(clu_dat_GO_res)){
  write.table(clu_dat_GO_res[[i]],paste0("~/Dropbox (PopGen)/Wei-Yun/project_aging/maSigPro_GO_v2/",names(clu_dat_GO_res)[i],"_GO_res.txt"),quote = F,
              row.names = F,sep="\t")
}
####maSigPro_cluster Rcistartget analysis####
# download.file("https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc_v10_clust/gene_based/dm6_v10_clust.genes_vs_motifs.rankings.feather",
#               destfile = "./dm6_v10_clust.genes_vs_motifs.rankings.feather")
ensembl=useMart("ensembl",dataset="dmelanogaster_gene_ensembl")
conv_background=ID_converter(ID=rownames(y),db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1]

motifranking=importRankings("./dm6-5kb-upstream-full-tx-11species.mc8nr.feather")
data(motifAnnotations_dmel_v8)
all_TFs=intersect(unique(motifAnnotations_dmel_v8$TF),conv_background)
all_TFs_tmp=ID_converter(all_TFs,db = ensembl,attributes = c("flybase_gene_id","external_gene_name","description"),
                         filters = "external_gene_name")
all_TFs_tmp=all_TFs_tmp[all_TFs_tmp$external_gene_name%in%all_TFs,]
# gene_summary_filtered_tmp=gene_summary_filtered[!duplicated(gene_summary_filtered$FB_ID_new),]
# rownames(gene_summary_filtered_tmp)=gene_summary_filtered_tmp$FB_ID_new
# all_TFs_out=cbind(all_TFs_tmp[,1:2],gene_summary_filtered_tmp[all_TFs_tmp$flybase_gene_id,1:4],
#                   sex_bias=apply(sapply(sex_bias_ID_new[c(1,2)],function(x) ifelse(all_TFs_tmp$flybase_gene_id%in%x,1,0)),1,function(x) ifelse(sum(x)==0,NA,names(which(x==1)))),
#                   evolution=apply(sapply(query_ID_new[-c(5,8,9,12)],function(x) ifelse(all_TFs_tmp$flybase_gene_id%in%x,1,0)),1,function(x) ifelse(sum(x)==0,NA,names(which(x==1)))),
#                   description=all_TFs_tmp[,3])

clu_ID[[5]]=clu_ID[[5]][-1]
clu_ID_conv=parallel::mclapply(clu_ID,function(x) ID_converter(ID=x,db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1],mc.cores = 12)
motifEnrichmentTable=parallel::mclapply(clu_ID_conv,function(x) {
  cisTarget(x,motifRankings = motifranking,motifAnnot = motifAnnotations_dmel_v8,nesThreshold = 3,nCores = 2,
            geneErnMethod = "iCisTarget",aucMaxRank = 0.01*ncol(motifranking))
},mc.cores = 12)
TFs=lapply(motifEnrichmentTable,function(x) {
  tmp=c(x$TF_highConf)
  genes <- gsub(" \\(.*\\). ", "; ", tmp, fixed=FALSE)
  genesSplit <- unique(unlist(strsplit(genes, "; ")))
  genesSplit <- unique(strsplit2(genesSplit, " ")[,1])
  return(genesSplit)
})


for (i in 1:12){
  tmp=which(sapply(clu_ID_conv,function(x) any(x%in%TFs[[i]])))
  print(c(i,tmp))
}

TFs_selected=lapply(TFs,function(x) x[x%in%unlist(clu_ID_conv)])

lapply(TFs_selected,function(x) sapply(x,function(a) which(sapply(clu_ID_conv,function(b) any(b%in%a)))))

TFoutTab=data.frame(TF=unique(unlist(TFs_selected)),
                    target=sapply(unique(unlist(TFs_selected)),function(x) paste(which(sapply(TFs_selected,function(a) any(a%in%x))),collapse = ",")),
                    pattern=sapply(unique(unlist(TFs_selected)),function(x) which(sapply(clu_ID_conv,function(a) any(a%in%x)))))
write.table(TFoutTab,"./TF_table.txt",quote = F,row.names = F)

dir.create("./maSigPro_GO_v2/")
for (i in 1:length(clu_dat_GO_res)){
  write.table(clu_dat_GO_res[[i]],paste0("./maSigPro_GO_v2/",names(clu_dat_GO_res)[i],"_GO_res.txt"),quote = F,
              row.names = F,sep="\t")
}
####maSigPro_tissue_enrichment####
flyatlas=read.table("/Volumes/Temp1/shengkai/common_info/fly_atlas_enrichment.table",header = T,stringsAsFactors = F)
flyatlas2=read.table("/Volumes/Temp1/shengkai/common_info/flyatlas2/version_18.05.25/flyaltlas2_log2fc.txt",header = T,stringsAsFactors = F)
tissues=read.delim("/Volumes/Temp1/shengkai/common_info/flyatlas2/version_18.05.25/Tissue.txt",header = T,stringsAsFactors = F)

p.val=c()
odds=c()
exp.num=c()
tissue_enriched_ID=list()
for(i in seq(4,36,2)){
  tissue_specific_ID=flyatlas$FB[flyatlas[,i]>2]
  p=c()
  o=c()
  e=c()
  g=list()
  for (j in 1:length(unique(cut))){
    p=c(p,fisher.test(cont_table(clu_ID[[j]],rownames(clusterdata),tissue_specific_ID),alternative = "greater")$p.value)
    o=c(o,fisher.test(cont_table(clu_ID[[j]],rownames(clusterdata),tissue_specific_ID),alternative = "greater")$estimate)
    g[[j]]=intersect(rownames(clusterdata)[idx],tissue_specific_ID)
    e=c(e,length(rownames(clusterdata)[idx])*sum(tissue_specific_ID%in%rownames(clusterdata))/length(rownames(clusterdata)))
  }  
  p.val=rbind(p.val,p)
  odds=rbind(odds,o)
  tissue_enriched_ID[[(i-2)/2]]=g
  exp.num=rbind(exp.num,e)
}
odds=odds[-c(7:9,13:14),]
exp.num=exp.num[-c(7:9,13:14),]
exp.num=round(exp.num,0)
obs.num= t(sapply(tissue_enriched_ID,function(x) sapply(x,length)))[-c(7:9,13:14),]
padj=matrix(p.adjust(p.val[-c(7:9,13:14),],method = "BH"),12,12)
names(tissue_enriched_ID)=strsplit2(colnames(flyatlas)[seq(4,36,2)],split = "_")[,1]
rownames(padj)=c("Br","Hd","Cr","Mg","Hg","Tb","Tg","Cs","Sg","Fb","Ey","Hr")
colnames(padj)=paste(1:12,table(cut))
rownames(odds)=c("Br","Hd","Cr","Mg","Hg","Tb","Tg","Cs","Sg","Fb","Ey","Hr")
colnames(odds)=paste(1:12,table(cut))
rownames(exp.num)=c("Br","Hd","Cr","Mg","Hg","Tb","Tg","Cs","Sg","Fb","Ey","Hr")
rownames(obs.num)=c("Br","Hd","Cr","Mg","Hg","Tb","Tg","Cs","Sg","Fb","Ey","Hr")

breaks=-log10(c(1,0.05,0.01,0.001,1e-1000))
color=c("lightyellow","gold","orange","firebrick")
pheatmap(-log10(padj),cluster_cols = F,cluster_rows = F,breaks = breaks,color = color,legend = F)

##watchout_don't run this line
#pheatmap(t(odds[,c(8,5,12,9)]),cluster_cols = F,cluster_rows = F,color = brewer.pal(5,"OrRd"),cex=1.5,display_numbers = T,
#         filename = "florida_result_fullset_v7/Figure1D.png",height = 8/2.54,width = 24/2.54,legend = F,
#         labels_row = c("F.up","F.down","M.up","M.down"))
#dev.off()

# write.csv(padj,"./florida_result_fullset_v7/enrichment_analysis/tissue_enrichment_padj.csv",quote = F)
# write.csv(odds,"./florida_result_fullset_v7/enrichment_analysis/tissue_enrichment_odds.csv",quote = F)
# write.csv(obs.num,"./florida_result_fullset_v7/enrichment_analysis/tissue_enrichment_obs.csv",quote = F)
# write.csv(exp.num,"./florida_result_fullset_v7/enrichment_analysis/tissue_enrichment_exp.csv",quote = F)


p.val2=c()
tissue_enriched_ID2=list()
for(i in 2:dim(flyatlas2)[2]){
  tissue_specific_ID2=flyatlas2$FBgnID[!is.na(flyatlas2[,i])&flyatlas2[,i]>1]
  p=c()
  g=list()
  for (j in 1:length(unique(cut))){
    # idx=cut==j
    p=c(p,fisher.test(cont_table(clu_ID[[j]],rownames(clusterdata),tissue_specific_ID2),alternative = "greater")$p.value)
    g[[j]]=intersect(clu_ID[[j]],tissue_specific_ID2)
  }  
  p.val2=rbind(p.val2,p)
  tissue_enriched_ID2[[i-1]]=g
}
padj2=matrix(p.adjust(p.val2[1:13,],method = "BH"),13,12)
names(tissue_enriched_ID2)=colnames(flyatlas2)[2:14]
rownames(padj2)=colnames(flyatlas2)[2:14]
colnames(padj2)=paste('Cluster',1:12,table(cut))
library(pheatmap)
breaks=-log10(c(1,0.05,0.01,0.001,1e-5,1e-10,1e-1000))
color=c("lightyellow","lightgoldenrod1","gold","orange","orangered","firebrick")
png("/Volumes/Temp2/shengkai/age_CGE/formal_experiment/maSigPro_tissue_enrichment.png",height = 16,width = 16,units = "cm",res = 600,pointsize = 10)
pheatmap(-log10(padj2),cluster_cols = F,cluster_rows = F,breaks = breaks,color = color,legend = F)
dev.off()

####quick check for Dagny####
dat.subset=count_dat_filtered[,substr(group,3,4)%in%c(3,5,7)]
group.subset=group[substr(group,3,4)%in%c(3,5,7)]
y.subset=DGEList(counts = dat.subset,group = group.subset)
y.subset=calcNormFactors(y.subset)

evo.s=strsplit2(group.subset,"_")[,1]
age.s=strsplit2(group.subset,"_")[,2]
ModelDesign.s=model.matrix(~evo.s+age.s+evo.s:age.s)
DGE.s=estimateDisp(y.subset,design = ModelDesign.s,robust = T)
GLM.s=glmFit(DGE.s,design = ModelDesign.s)
LRT_res_evo.s=glmLRT(GLM.s,coef = 2)
res_table_evo.s=LRT_res_evo.s$table
res_table_evo.s$padj=p.adjust(res_table_evo.s$PValue,method = "BH")
sum(res_table_evo.s$padj<0.05)

LRT_res_age.s=glmLRT(GLM.s,coef = 3:4)
res_table_age.s=LRT_res_age.s$table
res_table_age.s$padj=p.adjust(res_table_age.s$PValue,method = "BH")
sum(res_table_age.s$padj<0.05)

LRT_res_inter.s=glmLRT(GLM.s,coef = 5:6)
res_table_inter.s=LRT_res_inter.s$table
res_table_inter.s$padj=p.adjust(res_table_inter.s$PValue,method = "BH")
sum(res_table_inter.s$padj<0.05)

query_ID.s=list("evo_up"=rownames(y.subset)[res_table_evo.s$padj<0.05&res_table_evo.s$logFC>0],
                "evo_dn"=rownames(y.subset)[res_table_evo.s$padj<0.05&res_table_evo.s$logFC<0],
                "age_up"=rownames(y.subset)[res_table_age.s$padj<0.05&res_table_age.s$logFC.age.s5>0&res_table_age.s$logFC.age.s7>0],
                "age_dn"=rownames(y.subset)[res_table_age.s$padj<0.05&res_table_age.s$logFC.age.s5<0&res_table_age.s$logFC.age.s7<0],
                "inter"=rownames(y.subset)[res_table_inter.s$padj<0.05])
p.val2.s=c()
for(i in 2:dim(flyatlas2)[2]){
  tissue_specific_ID2=flyatlas2$FBgnID[!is.na(flyatlas2[,i])&flyatlas2[,i]>1]
  p=c()
  for (j in 1:length(query_ID.s)){
    p=c(p,fisher.test(cont_table(query_ID.s[[j]],rownames(y.subset),tissue_specific_ID2),alternative = "greater")$p.value)
  }  
  p.val2.s=rbind(p.val2.s,p)
}  
padj2.s=matrix(p.adjust(p.val2.s[1:13,],method = "BH"),13,5)
rownames(padj2.s)=colnames(flyatlas2)[2:14]
colnames(padj2.s)=paste(names(query_ID.s),sapply(query_ID.s,length))
breaks=-log10(c(1,0.05,0.01,0.001,1e-5,1e-10,1e-1000))
color=c("lightyellow","lightgoldenrod1","gold","orange","orangered","firebrick")
pheatmap(-log10(padj2.s),cluster_cols = F,cluster_rows = F,breaks = breaks,color = color,legend = F)

mg_expressed_gene=flyatlas2$FBgnID[!is.na(flyatlas2[,7])&flyatlas2[,7]>1]
Reduce(intersect,list(query_ID.s$evo_dn,query_ID.s$age_dn,mg_expressed_gene))

expr.dat.s=log(cpm(y.subset))
scl.expr.dat.s=t(apply(expr.dat.s,1,function(x) (x-mean(x))/sd(x)))
avg.expr.dat.s=apply(scl.expr.dat.s,1,function(x) tapply(x,group.subset,mean))
cand=sample(Reduce(intersect,list(query_ID.s$evo_dn,query_ID.s$age_dn,mg_expressed_gene)),10)

plot(NA,xlim=c(1,6),ylim=c(-2,2),xlab='age (day)',ylab='scaled expression',xaxt='n')
axis(1,at=1:6,labels = rep(seq(3,7,2),2))
for (i in 1:length(cand)){
  lines(1:6,t(avg.expr.dat.s)[cand[i],],lty=2,lwd=0.75)
  lines(1:3,t(avg.expr.dat.s)[cand[i],1:3],col=alpha("forestgreen",0.5),lwd=1.5)
  lines(4:6,t(avg.expr.dat.s)[cand[i],4:6],col=alpha("salmon",0.5),lwd=1.5)
}

plot(res_table_age.s$logFC.age.s7,res_table_evo.s$logFC,pch=1,asp=1,xlim=c(-1,1),ylim=c(-1,1),col="grey",cex=0.5)
abline(h=0,col='blue',lty=3)
abline(v=0,col='blue',lty=3)
points(res_table_age.s[mg_expressed_gene,]$logFC.age.s7,res_table_evo.s[mg_expressed_gene,]$logFC,pch=19,cex=0.5)
points(res_table_age.s[cand,]$logFC.age.s7,res_table_evo.s[cand,]$logFC,pch=19,cex=0.5,col='red')

plotMDS(y.subset[rownames(y.subset)%in%mg_expressed_gene,])
pca.s=prcomp(t(log(cpm(y.subset))))
pca.s.mg=prcomp(t(log(cpm(y.subset)[rownames(y.subset)%in%mg_expressed_gene,])))
png("/Volumes/Temp2/shengkai/age_CGE/formal_experiment/mg_expressed_gene_PCA_Dagny.png",
    width = 17.4,height = 8.7,units = "cm",res=600,pointsize = 8)
par(mfrow=c(1,2),mar=c(4,4,2,2))
plot(pca.s$x,pch=ifelse(evo.s%in%'B',1,19),col=as.factor(age.s),asp=1,ylim=c(-12,12))
legend("topleft",pch = c(1,19),legend = c("Anc.", "Evolved"),bty='n')
plot(pca.s.mg$x,pch=ifelse(evo.s%in%'B',1,19),col=as.factor(age.s),asp=1,ylim=c(-6,6))
legend("bottomright",fill = c(1,2,3),legend = c(3,5,7),bty='n')
dev.off()

####old stuffs####
#linear modeling and DE analysis
ModelDesign=model.matrix(~0+group)
DGE=estimateDisp(y,design = ModelDesign,robust = T)
GLM=glmFit(DGE,design = ModelDesign)
my.contrasts=makeContrasts(ageB=groupB_50-groupB_05,
                           ageH=groupH_50-groupH_05,
                           evo05=groupH_05-groupB_05,
                           evo50=groupH_50-groupB_50,
                           inter=(groupH_50-groupB_50)-(groupH_05-groupB_05),
                           levels=ModelDesign)

LRT_list=list()
for (i in colnames(my.contrasts)){
  LRT_list[[i]]=glmLRT(GLM,contrast = my.contrasts[,i])
}
res_list=lapply(LRT_list,function(x) x$table)
res_list=lapply(res_list,function(x) 
  {x$padj=p.adjust(x$PValue,method = "BH")
  return(x)})
sapply(res_list,function(x) sum(x$padj<0.05))

png("./scatter_plot_trial.png",width = 12,height = 12, units = "cm",res = 600,pointsize = 8)
plot(res_list$evo05$logFC,res_list$evo50$logFC,xlab="evolution (5-day-old)",ylab="evolution (50-day-old)",
     pch=19,asp=1)
abline(h=0,v=0,col="grey70",lty=2,lwd=2)
abline(a=0,b=1,col="grey70",lty=2,lwd=2)
abline(a=0,b=-1,col="grey70",lty=2,lwd=2)
points(res_list$evo05$logFC[res_list$inter$padj<0.05],res_list$evo50$logFC[res_list$inter$padj<0.05],col="red",pch=19)
text(2,-2,labels = expression(rho==0.12))
dev.off()
cor.test(res_list$evo05$logFC,res_list$evo50$logFC,method="spearman")

png("./scatter_plot2_trial.png",width = 12,height = 12, units = "cm",res = 600,pointsize = 8)
plot(res_list$ageB$logFC,res_list$ageH$logFC,xlab="Aging in anc. population",ylab="Aging in hot evolved population",
     pch=19,asp=1)
abline(h=0,v=0,col="grey70",lty=2,lwd=2)
abline(a=0,b=1,col="grey70",lty=2,lwd=2)
abline(a=0,b=-1,col="grey70",lty=2,lwd=2)
points(res_list$ageB$logFC[res_list$inter$padj<0.05],res_list$ageH$logFC[res_list$inter$padj<0.05],col="red",pch=19)
text(5,-2,labels = expression(rho==0.87))
dev.off()
cor.test(res_list$ageB$logFC,res_list$ageH$logFC,method="spearman")

background=rownames(y)
query_ID=list(evo05_up=background[res_list$evo05$padj<0.05&res_list$evo05$logFC>0],
              evo05_dn=background[res_list$evo05$padj<0.05&res_list$evo05$logFC<0],
              evo50_up=background[res_list$evo50$padj<0.05&res_list$evo50$logFC>0],
              evo50_dn=background[res_list$evo50$padj<0.05&res_list$evo50$logFC<0],
              ageB_up=background[res_list$ageB$padj<0.05&res_list$ageB$logFC>0],
              ageB_dn=background[res_list$ageB$padj<0.05&res_list$ageB$logFC<0],
              ageH_up=background[res_list$ageH$padj<0.05&res_list$ageH$logFC>0],
              ageH_dn=background[res_list$ageH$padj<0.05&res_list$ageH$logFC<0],
              inter_50=background[res_list$inter$padj<0.05&res_list$inter$logFC>0],
              inter_05=background[res_list$inter$padj<0.05&res_list$inter$logFC<0],
              cons_evo_up=background[res_list$evo05$padj<0.05&res_list$evo05$logFC>0&res_list$evo50$padj<0.05&res_list$evo50$logFC>0],
              cons_evo_dn=background[res_list$evo05$padj<0.05&res_list$evo05$logFC<0&res_list$evo50$padj<0.05&res_list$evo50$logFC<0],
              anta_evo_50=background[res_list$evo05$padj<0.05&res_list$evo05$logFC<0&res_list$evo50$padj<0.05&res_list$evo50$logFC>0],
              anta_evo_05=background[res_list$evo05$padj<0.05&res_list$evo05$logFC>0&res_list$evo50$padj<0.05&res_list$evo50$logFC<0],
              evo05_s_up=background[res_list$evo05$padj<0.05&res_list$evo05$logFC>0&res_list$evo50$padj>0.05],
              evo05_s_dn=background[res_list$evo05$padj<0.05&res_list$evo05$logFC<0&res_list$evo50$padj>0.05],
              evo50_s_up=background[res_list$evo50$padj<0.05&res_list$evo50$logFC>0&res_list$evo05$padj>0.05],
              evo50_s_dn=background[res_list$evo50$padj<0.05&res_list$evo50$logFC<0&res_list$evo05$padj>0.05])
sapply(query_ID,length)

sapply(res_list,function(x) x["FBgn0014859",])#Hr38
sapply(res_list,function(x) x["FBgn0005626",])#ple
sapply(res_list,function(x) x["FBgn0000422",])#Ddc
sapply(res_list,function(x) x["FBgn0019643",])#Dat
sapply(res_list,function(x) x["FBgn0260964",])#Vmat
sapply(res_list,function(x) x["FBgn0034136",])#DAT
sapply(res_list,function(x) x["FBgn0011656",])#mef2

#enrichment analysis
flyatlas=read.table("/Volumes/Temp1/shengkai/common_info/fly_atlas_enrichment.table",header = T,stringsAsFactors = F)
flyatlas2=read.table("/Volumes/Temp1/shengkai/common_info/flyatlas2/version_18.05.25/flyaltlas2_log2fc.txt",header = T,stringsAsFactors = F)
tissues=read.delim("/Volumes/Temp1/shengkai/common_info/flyatlas2/version_18.05.25/Tissue.txt",header = T,stringsAsFactors = F)

ensembl=useMart("ensembl",dataset="dmelanogaster_gene_ensembl")
background_new=background
latest_FB_ID=read.delim("/Volumes/Temp1/shengkai/common_info/fbgn_annotation_ID_fb_2018_01_R.txt",header = T,stringsAsFactors = F,sep="\t")
for(i in 1:dim(latest_FB_ID)[1]){
  background_new[which(background%in%strsplit2(latest_FB_ID[i,4],","))]=latest_FB_ID[i,3]
}
query_ID_new=lapply(query_ID,function(x) {
  for(i in 1:dim(latest_FB_ID)[1]){
    x[which(x%in%strsplit2(latest_FB_ID[i,4],","))]=latest_FB_ID[i,3]
  }
  return(x)
})
conv_query=lapply(query_ID_new,function(x) ID_converter(ID = x,db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1])
conv_background=ID_converter(ID=background_new,db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1]

for (i in names(query_ID_new)){
  write.table(query_ID_new[[i]],paste0("./GOI/GOI_",i,".txt"),sep="\t",quote = F,row.names = F)
}

for (i in names(query_ID_new)){
  write.table(conv_query[[i]],paste0("./GOI/converted/",i,".txt"),sep="\t",quote = F,row.names = F)
}

#GO
GO_res_table=list()
for (i in names(query_ID_new)){
  tmp=factor(as.integer(background_new%in%query_ID_new[[i]]))
  names(tmp)=background_new
  tgd=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
  resTopGO.classic=runTest(tgd, algorithm = "classic", statistic = "Fisher")
  resTopGO.weight01=runTest(tgd, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)
  #resTopGO.ParentChild=runTest(tgd, algorithm = "parentchild", statistic = "Fisher")
  #tmp_res=GenTable(tgd,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,Fisher.parentchild=resTopGO.ParentChild,orderBy = "Fisher.weight01",ranksOf="Fisher.parentchild",topNodes=length(resTopGO.classic@score),numChar=100)
  GO_res_table[[i]]=tmp_res
}
GO_res_table=lapply(GO_res_table,function(x) {
  x$Fisher.weight01[x$Fisher.weight01=="< 1e-30"]=1e-30
  return(x)})

for (i in names(query_ID_new)){
  write.table(GO_res_table[[i]],paste0("./GO_enrichment/GO_",i,".txt"),sep="\t",quote = F,row.names = F)
}

p.val2=c()
tissue_enriched_ID2=list()
for(i in 2:dim(flyatlas2)[2]){
  tissue_specific_ID2=flyatlas2$FBgnID[!is.na(flyatlas2[,i])&flyatlas2[,i]>1]
  p=c()
  g=list()
  for (j in 1:length(query_ID_new)){
    p=c(p,fisher.test(cont_table(query_ID_new[[j]],background,tissue_specific_ID2),alternative = "greater")$p.value)
    g[[j]]=intersect(query_ID_new[[j]],tissue_specific_ID2)
  }  
  p.val2=rbind(p.val2,p)
  names(g)=names(query_ID_new)
  tissue_enriched_ID2[[i-1]]=g
}
padj2=matrix(p.adjust(p.val2[1:13,],method = "BH"),13,18)
names(tissue_enriched_ID2)=colnames(flyatlas2)[2:37]
rownames(padj2)=colnames(flyatlas2)[2:14]
colnames(padj2)=paste(names(query_ID),sapply(query_ID,length))
library(pheatmap)
breaks=-log10(c(1,0.05,0.01,0.001,1e-5,1e-10,1e-1000))
color=c("lightyellow","lightgoldenrod1","gold","orange","orangered","firebrick")
png("/Volumes/Temp2/shengkai/age_CGE/heatmap_tissue_enrichment_padj2.png",height = 16,width = 16,units = "cm",res = 600,pointsize = 10)
pheatmap(-log10(padj2),cluster_cols = F,cluster_rows = F,breaks = breaks,color = color,legend = F)
dev.off()

GOI_ACP=read.table("/Volumes/Temp1/shengkai/CHC/ACP_flybase.txt",stringsAsFactors = F)[,1]
GOI_SFP=read.table("/Volumes/Temp1/shengkai/CHC/SFP_flybase.txt",stringsAsFactors = F)[,1]
fisher.test(cont_table(query_ID_new$inter_05,background,GOI_ACP),alternative = "greater")
fisher.test(cont_table(query_ID_new$inter_05,background,GOI_SFP),alternative = "greater")

png("ACP_boxplot_trial.png",width = 16,height = 12,units = "cm",res = 600,pointsize = 8)
boxplot(t(apply(apply(log10(cpm(y)[intersect(GOI_ACP,background),]),1,function(x) tapply(x,group,mean)),2,scale)),
        names=c("Base (5-d-old)","Base (50-d-old)","Hot evolved (5-d-old)","Hot evolved (50-d-old)"),
        xlab="Population",ylab="Normalized expression")
dev.off()

GOI_reproduction=intersect(genesInTerm(tgd,"GO:0032504")[[1]],query_ID_new$inter_05)
png("reproductive_genes_boxplot_trial.png",width = 16,height = 12,units = "cm",res = 600,pointsize = 8)
boxplot(t(apply(apply(log10(cpm(y)[intersect(GOI_reproduction,background),]),1,function(x) tapply(x,group,mean)),2,scale)),
        names=c("Base (5-d-old)","Base (50-d-old)","Hot evolved (5-d-old)","Hot evolved (50-d-old)"),
        xlab="Population",ylab="Normalized mean expression")
dev.off()

dop_GOI=rev(c("FBgn0014859","FBgn0005626","FBgn0000422","FBgn0019643","FBgn0260964","FBgn0034136",
              "FBgn0015129","FBgn0053517","FBgn0035538"))
dop_GOI_s=rev(c("Hr38","Ple","Ddc","Dat","Vmat","DAT","Dop1R2","Dop2R","DopEcR"))

png("DPN_boxplot_trial.png",width = 16,height = 12,units = "cm",res = 600,pointsize = 8)
boxplot(t(apply(apply(log10(cpm(y)[intersect(dop_GOI,background),]),1,function(x) tapply(x,group,mean)),2,scale)),
        names=c("Base (5-d-old)","Base (50-d-old)","Hot evolved (5-d-old)","Hot evolved (50-d-old)"),
        xlab="Population",ylab="Normalized expression")
dev.off()


GOI_circadian=intersect(genesInTerm(tgd,"GO:0045475")[[1]],query_ID_new$inter_50)
png("Circadian_boxplot_trial.png",width = 16,height = 12,units = "cm",res = 600,pointsize = 8)
boxplot(t(apply(apply(log10(cpm(y)[intersect(GOI_circadian,background),]),1,function(x) tapply(x,group,mean)),2,scale)),
        names=c("Base (5-d-old)","Base (50-d-old)","Hot evolved (5-d-old)","Hot evolved (50-d-old)"),
        xlab="Population",ylab="Normalized expression")
dev.off()
