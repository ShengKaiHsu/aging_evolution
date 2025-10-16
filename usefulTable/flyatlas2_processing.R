#flyatlas2 DB processing
##SKHsu_20180227
##SKHsu_20180702_modified (v18.5.25)

rm(list=ls())
setwd("/Volumes/Temp1/shengkai/common_info/flyatlas2/version_18.05.25/")
geneFPKM=read.table("./geneFPKM.txt",header=T,stringsAsFactors = F)
trans_tissue=read.delim("./Tissue.txt",header=T,stringsAsFactors = F)

meta_table=c()
for (i in 1:length(unique(geneFPKM$FBgn))){
  print(i)
  meta_table=rbind(meta_table,c(unique(geneFPKM$FBgn)[i],geneFPKM[geneFPKM$FBgn%in%unique(geneFPKM$FBgn)[i],]$FPKM))
}
meta_table=data.frame(meta_table[,1],apply(meta_table[,-1],2,function(x) as.numeric(x)))
colnames(meta_table)=c("FBgnID",with(trans_tissue,paste(Stage,Sex,Abbreviation,sep="_")[trans_tissue$TissueID%in%unique(geneFPKM$TissueID)]))

meta_table_sd=c()
for (i in 1:length(unique(geneFPKM$FBgn))){
  print(i)
  meta_table_sd=rbind(meta_table_sd,c(unique(geneFPKM$FBgn)[i],geneFPKM[geneFPKM$FBgn%in%unique(geneFPKM$FBgn)[i],]$SD))
}
meta_table_sd=data.frame(meta_table_sd[,1],apply(meta_table_sd[,-1],2,function(x) as.numeric(x)))
colnames(meta_table_sd)=c("FBgnID",with(trans_tissue,paste(Stage,Sex,Abbreviation,sep="_")[trans_tissue$TissueID%in%unique(geneFPKM$TissueID)]))

meta_table_fc=data.frame(FBgnID=meta_table$FBgnID)
for (i in 1:3){
  tmp=meta_table[,grep(unique(substr(colnames(meta_table[,-1]),1,7))[i],colnames(meta_table))]
  fc.tmp=log2(t(apply(tmp,1,function(x) x/x[length(x)])))
  meta_table_fc=data.frame(meta_table_fc,fc.tmp)
}
write.table(meta_table,"./flyaltlas2_avg_expression.txt",row.names = F,quote = F)
write.table(meta_table_fc,"./flyaltlas2_log2fc.txt",row.names = F,quote = F)
write.table(meta_table_sd,"./flyaltlas2_expression_sd.txt",row.names = F,quote = F)


