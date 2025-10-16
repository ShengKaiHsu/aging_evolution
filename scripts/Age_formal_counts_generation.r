#count_data_generation_complete_dataset
#SKHsu_2017_11_14
#modified from Ana marija´s code

rm(list=ls())
gc()
library("edgeR")
library("gplots")
library("RColorBrewer")
library(ggplot2)
library(matrixStats)
library(limma)

directory_path="/Volumes/Temp2/shengkai/age_CGE/formal_experiment/ageing_formal_output/Rsubread/"
setwd(directory_path)
all_files=list.files(path=directory_path,pattern=".tsv")
read_all=lapply(all_files,read.csv,sep="\t",header=TRUE)#,row.names=1)

temp=as.data.frame(read_all[1]);temp$gene_id=as.character(row.names(temp));temp=temp[,c(2,1)]
names(temp)=c("gene_id",all_files[1])
counts=temp

dim(counts)
names(counts)
head(counts)

for(i in 2:length(read_all)){
  
  temp=as.data.frame(read_all[i]);temp$gene_id=as.character(row.names(temp));temp=temp[,c(2,1)]
  names(temp)=c("gene_id",all_files[i])
  temp$gene_id=as.character(temp$gene_id)
  
  counts=merge(counts,temp,all.x=T,all.y=T)
  print(dim(counts))
}

row.names(counts)=counts$gene_id
counts=subset(counts,!is.na(counts[,2]))
counts=counts[,2:dim(counts)[2]]
dim(counts)


ntrans=read.table("/Volumes/Temp2/shengkai/age_CGE/formal_experiment/Rsubread_translate_Age_formal.txt")
original=substr(names(counts),1,12)
modified=original
for(i in 1:length(original)) modified[i]=as.character(ntrans$V2)[which(ntrans$V1==original[i])[1]]

modified
names(counts)=modified
counts_dsim=counts[,!is.na(colnames(counts))]
write.csv(counts_dsim,"/Volumes/Temp2/shengkai/age_CGE/formal_experiment/age_formal_fullset_readcounts.csv",quote=F)

#hot=counts_dsim[,grep("2818", colnames(counts_dsim))]
#hot=hot[,order(colnames(hot)),]
#hotc=cpm(hot, normalized.lib.sizes=T)
#hotc=as.data.frame(hotc)
#hotexx=hot[rownames(hotc[order(apply(hotc,1,mean), decreasing=T),][1:round(length(hotc[,1])*0.7),]),]


