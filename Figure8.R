 
################################################################################
# Figure 8B DEGs
################################################################################
# MAST was used to calculate sex-related differential genes in ECs from different tissues
library(MAST);library(Seurat);library(dplyr)
#load("/home/niuruize2/aging/PFC/merge/scRNA_PFC_seurat.RData")
Idents(scRNA)="Sex"
table(Idents(scRNA))
scRNA1 = scRNA
fData = data.frame(symbolid=rownames(scRNA1),primerid=rownames(scRNA1))
rownames(fData)=fData$symbolid
cData = scRNA1@meta.data
cData$wellKey <- rownames(cData)
sca = FromMatrix(as.matrix(scRNA1@assays$RNA@data), cData = cData,fData = fData)
rm(scRNA1)
gc()
dim(sca)
table(colData(sca)$Sex)
cond<-factor(colData(sca)$Sex)
cond<-relevel(cond,"M")
colData(sca)$condition<-cond

colData(sca)$Sex=factor(colData(sca)$Sex)
colData(sca)$Age=as.numeric(colData(sca)$Age)
colData(sca)$percent.mt=as.numeric(colData(sca)$percent.mt)
colData(sca)$percent.rb=as.numeric(colData(sca)$percent.rb)
colData(sca)$nCount_RNA=as.numeric(colData(sca)$nCount_RNA)
#colData(sca)$datasets=factor(colData(sca)$datasets)

## (1) 校正cngeneson协变量:默认参数

zlmCond <- zlm(~condition + nCount_RNA + percent.mt + percent.rb + Age, sca, method="bayesglm", ebayes=TRUE)
summaryCond <- summary(zlmCond,doLRT='conditionF')
summaryDt <- summaryCond$datatable
levels(summaryDt$contrast)

#整理结果
df_pval = summaryDt %>% 
  dplyr::filter(contrast=='conditionF') %>% 
  dplyr::filter(component=='H') %>% 
  dplyr::select(primerid, `Pr(>Chisq)`)

df_logfc = summaryDt %>% 
  dplyr::filter(contrast=='conditionF') %>% 
  dplyr::filter(component=='logFC') %>% 
  dplyr::select(primerid, coef, ci.hi, ci.lo)

df_stat = dplyr::inner_join(df_logfc, df_pval) %>% 
  dplyr::rename("symbol"="primerid") %>% 
  dplyr::rename("pval"="Pr(>Chisq)","logFC"="coef") %>% 
  dplyr::mutate("fdr" = p.adjust(pval)) %>% 
  dplyr::arrange(fdr)
head(df_stat)

df_stat$FC<-10^(abs(df_stat$logFC))
df_stat$FC<-ifelse(df_stat$logFC>0,df_stat$FC*(1),df_stat$FC*-1) 
df_stat$log2FC <- log2(abs(df_stat$FC))
df_stat$log2FC <- ifelse(df_stat$FC>0,df_stat$log2FC*(1),df_stat$log2FC*-1)
write.csv(df_stat,"Brain_F_M_DEGs.csv")

################################################################################
# Figure 8C Overlap of genes associated with sex from different tissues
################################################################################
library(dplyr);library(GeneOverlap);data(GeneOverlap)

#Select variables
variable_list=list("Group")

for (k in 1:length(variable_list)) {
  variable=variable_list[[k]]
  
  setwd("/Users/niuruize/Downloads/scRNA/EC/EC/10_sex")
  
  celltype = c("Aorta","Bladder","Brain","Heart","Intestinal","Kidney","Liver","Lung","Muscle","Pancreas","PIT","Retina","SAT","SC","Skin","VAT")
  
  #List significantly upregulated genes
  
  Marker_genes_list<-list()
  for (i in 1:length(celltype)) {
    ct1=celltype[i]
    ct1=gsub("/","_",ct1)
    input_file=paste0(ct1,"_F_M_DEGs",".csv")
    topgenes=read.csv(input_file,header=TRUE,row.names=1)
    Marker_genes_full=topgenes
    Marker_genes_full_significant=subset(Marker_genes_full,fdr < 0.05)
    Marker_genes_full_significant=subset(Marker_genes_full_significant,logFC > 0)
    Marker_genes=Marker_genes_full_significant[,1:2]
    Marker_genes[] <- lapply(Marker_genes, as.character)
    if (nrow(Marker_genes)>0){
      Marker_genes_list[[i]]<-Marker_genes[[1]]
      names(Marker_genes_list)[i]<- celltype[i]
    }
  }
  
  #Remove empty elements
  Marker_genes_list=Marker_genes_list[lengths(Marker_genes_list) != 0]
  
  for (c in 1:length(celltype)) {
    ct1=celltype[[c]]
    ct1_name=gsub("/","_",ct1)
    for (b in 1:length(celltype)){
      ct2=celltype[[b]]
      ct2_name=gsub("/","_",ct2)
      
      input_file1=paste0(ct1_name,"_F_M_DEGs",".csv")
      input_file2=paste0(ct2_name,"_F_M_DEGs",".csv")
      
      topgenes1=read.csv(input_file1,header=TRUE,row.names=1)
      
      topgenes2=read.csv(input_file2,header=TRUE,row.names=1)
      
      All_genes=full_join(topgenes1,topgenes2,by="symbol")
      
      genome_size=nrow(All_genes)
      
      gom.obj <- newGOM(Marker_genes_list,Marker_genes_list, genome.size = genome_size)
      GeneOverlap_Matrix=getMatrix(gom.obj, name="pval")
      GeneOverlap_Matrix_df=as.data.frame(GeneOverlap_Matrix)
      
      myData=GeneOverlap_Matrix_df[ct1,ct2,drop=FALSE]
      
      if (b==1){
        results_int=myData
      }
      if(b>1){
        results_int=cbind(results_int,myData)
      }
    }
    if (c==1){
      results_table=results_int
    }
    if (c>1){
      results_table=rbind(results_table,results_int)
    }
  }
  
  n_rows=nrow(results_table)
  n_columns=ncol(results_table)
  factor=n_rows*n_columns
  pvalue_adjusted_matrix <- apply(results_table,2,p.adjust,method="bonferroni", n = factor)
  pvalue_adjusted_matrix_minusLog10=-log10(pvalue_adjusted_matrix)
  pvalue_adjusted_matrix_minusLog10[is.infinite(pvalue_adjusted_matrix_minusLog10)] <- 300
  
  setwd("/Users/niuruize/Downloads/scRNA/EC/EC/10_sex")
  output_file=paste0("EC_sex_up_P",".csv")
  write.csv(pvalue_adjusted_matrix_minusLog10,file=output_file)
  
  
  GeneOverlap_Matrix=getMatrix(gom.obj, name="odds.ratio")
  GeneOverlap_Matrix[is.infinite(GeneOverlap_Matrix)] <- 20
  odds_ratio_Matrix = as.data.frame(GeneOverlap_Matrix)
  output_file=paste0("EC_sex_up_OR",".csv")
  write.csv(odds_ratio_Matrix,file=output_file)
}

library(pheatmap)
library(ggplot2)
r2 <- read.csv('EC_sex_up_OR.csv',header=TRUE,row.names=1)
p2 <- read.csv('EC_sex_up_P.csv',header=TRUE,row.names=1)

r2[r2 >= 10] <- 10

p2[p2 <= 10] <- ""
p2[p2 >= 100] <- "**"
p2[p2 > 10] <- "*"

#mycol<-colorRampPalette(c( "#0f86a9", "white", "#ed8b10"))(200)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)

#绘制热图；
p1 = pheatmap(r2,scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=14,fontsize_row=8,fontsize_col=9,cellwidth=15,
              cellheight=15,cluster_rows=T,cluster_cols=T,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "EC_sex_overlap_up1.pdf", p1,device = 'pdf', width = 15, height = 15, units = 'cm')
p1 = pheatmap(r2,scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=14,fontsize_row=8,fontsize_col=9,cellwidth=15,
              cellheight=15,cluster_rows=F,cluster_cols=F,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "EC_sex_overlap_up2.pdf", p1,device = 'pdf', width = 15, height = 15, units = 'cm')

################################################################################
# Figure 8E GO analysis
################################################################################
library("clusterProfiler");library("org.Hs.eg.db");library("enrichplot");library("ggplot2");library("org.Hs.eg.db")
##
#--------------------------------------------------------------------Dowm---------------------------------------------------
celltype = c("Aorta","Bladder","Brain","Heart","Intestinal","Kidney","Liver","Lung","Muscle","Pancreas","Retina","SAT","SC","Skin","Testis","Thyroid","Uterus","VAT")
for (i in 1:length(celltype)) {
  ct1=celltype[i]
  ct1=gsub("/","_",ct1)
  input_file=paste0("/Users/niuruize/Downloads/scRNA/EC/EC/10_sex/",ct1,"_F_M_DEGs",".csv")
  topgenes = read.csv(input_file,row.names=1)
  topgenes <- subset(topgenes, topgenes$fdr < 0.05 & topgenes$log2FC < 0)
  df <- data.frame(gene=topgenes$symbol)
  genes=as.vector(df[,1])
  entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrezIDs <- as.character(entrezIDs)
  id=cbind(df,1,entrezID=entrezIDs)
  colnames(id)=c("symbol","logFC","entrezID")
  rt=id
  rt=rt[is.na(rt[,"entrezID"])==F,]
  gene=rt$entrezID
  kk <- enrichGO(gene = gene,ont = "BP",OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, readable =T)
  out_file=paste0(ct1,"_down_BP.csv")
  write.csv(kk@result,out_file, row.names = F) 
}
# merge
celltype = c("Aorta","Bladder","Brain","Heart","Intestinal","Kidney","Liver","Lung","Muscle","Pancreas","Retina","SAT","SC","Skin","Testis","Thyroid","Uterus","VAT")
# down
GO_down <- data.frame(GO = numeric(0))
for (i in 1:length(celltype)) {
  ct1=celltype[i]
  ct1=gsub("/","_",ct1)
  input_file=paste0(ct1,"_down_BP.csv")
  topGO = read.csv(input_file,row.names=1)
  topGO <- subset(topGO, topGO$p.adjust < 0.05)
  topGO <- data.frame(topGO$Description, -log10(topGO$p.adjust))
  names(topGO) <- c("GO", celltype[i])
  GO_down <- merge(GO_down, topGO,by="GO",all=TRUE)
}
write.csv(GO_down,"GO_down.csv", row.names = F) 

#--------------------------------------------------------------------Up---------------------------------------------------
##up
celltype = c("Aorta","Bladder","Brain","Heart","Intestinal","Kidney","Liver","Lung","Muscle","Pancreas","PIT","Retina","SAT","SC","Skin","Testis","Thyroid","Uterus","VAT")
celltype = c("Pancreas","PIT","Retina","SAT","SC","Skin","Testis","Thyroid","Uterus","VAT")

for (i in 1:length(celltype)) {
  ct1=celltype[i]
  ct1=gsub("/","_",ct1)
  input_file=paste0("/Users/niuruize/Downloads/scRNA/EC/EC/10_sex/",ct1,"_F_M_DEGs",".csv")
  topgenes = read.csv(input_file,row.names=1)
  topgenes <- subset(topgenes, topgenes$fdr < 0.05 & topgenes$log2FC > 0)
  df <- data.frame(gene=topgenes$symbol)
  genes=as.vector(df[,1])
  entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrezIDs <- as.character(entrezIDs)
  id=cbind(df,1,entrezID=entrezIDs)
  colnames(id)=c("symbol","logFC","entrezID")
  rt=id
  rt=rt[is.na(rt[,"entrezID"])==F,]
  gene=rt$entrezID
  kk <- enrichGO(gene = gene,ont = "BP",OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, readable =T)
  out_file=paste0(ct1,"_up_BP.csv")
  write.csv(kk@result,out_file, row.names = F) 
}

# merge
celltype = c("Aorta","Bladder","Brain","Heart","Intestinal","Kidney","Liver","Lung","Muscle","Pancreas","PIT","Retina","SAT","SC","Skin","Testis","Thyroid","Uterus","VAT")
##down
GO_up <- data.frame(GO = numeric(0))
for (i in 1:length(celltype)) {
  ct1=celltype[i]
  ct1=gsub("/","_",ct1)
  input_file=paste0(ct1,"_up_BP.csv")
  topGO = read.csv(input_file,row.names=1)
  topGO <- subset(topGO, topGO$p.adjust < 0.05)
  topGO <- data.frame(topGO$Description, -log10(topGO$p.adjust))
  names(topGO) <- c("GO", celltype[i])
  GO_up <- merge(GO_up, topGO,by="GO",all=TRUE)
}
write.csv(GO_up,"GO_up.csv", row.names = F) 
# plot the bar graph using Prism 9

################################################################################
# Figure 8F AutoML-H2O
################################################################################
library(h2o)
h2o.init()
library(Seurat);library(dplyr);library(randomForest);library(caret);library(pROC);library(caret);library(caret)
library(modeldata);library(tidymodels);library(recipes) 
# load data
load("Human_EC.RData")
# Extracting feature information
table(scRNA$Group)
scRNA = scRNA[,scRNA$Group %in% c("Adult","Aging")]
scRNA = scRNA[,scRNA$Sex %in% c("F")]
DefaultAssay(scRNA) <- "RNA"
sce.markers <- read.csv("markers.txt",sep = "")
Idents(scRNA) <- "Group"
#sce <- subset(scRNA, subset = nFeature_RNA > 500)
#sce <- subset(sce,downsample = 200)
sce <- scRNA
table(Idents(sce))
sce <- ScaleData(sce,features = unique(sce.markers$gene)) 
t_expr <- t(as.matrix(sce@assays$RNA@scale.data))
dim(t_expr) 
t_expr[1:4,1:4]

# train data and test data
inTrain <- createDataPartition(y= Idents(sce),p=0.25,list=F)
test_expr <- t_expr[inTrain,]
train_expr <- t_expr[-inTrain,]
test_y <- Idents(sce)[inTrain]
train_y <- Idents(sce)[-inTrain]
#save(test_y,train_y, test_expr,train_expr,file = 'input100.Rdata') 
train_expr[1:4,1:4]
table(train_y)
table(test_y)
table(scRNA$Group)

# Building data
test_expr = cbind(data.frame(test_expr),data.frame(Group=test_y))
train_expr = cbind(data.frame(train_expr),data.frame(Group=train_y))
#train_expr = cbind(data.frame(t_expr),data.frame(Group=Idents(sce)))
rec <- recipe(Group~.,train_expr) %>%           
  step_dummy(all_nominal_predictors()) %>%           
  prep() %>%           
  bake(new_data=NULL) %>%          
  as.h2o() 
rec$Group <- h2o.asfactor(rec$Group)
y <- "Group" 
x <- setdiff(names(rec),y) 

# Building model
am <- h2o.automl(x,y, training_frame = rec,max_models = 20, seed=1234)
save(am,file="h2o_automl.RData")
#am <- h2o.automl(y="Group", training_frame = train_expr, max_models = 10)

#
b <- h2o.get_leaderboard(am) 
b
best <- h2o.get_best_model(am)
#best <- h2o.get_best_model(am,"gbm")
perf <- h2o.performance(best)
perf
Confusion_Matrix <- as.matrix(perf@metrics$cm$table[1:2,1:2])
Confusion_Matrix <- prop.table(Confusion_Matrix, margin = 1)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap::pheatmap(Confusion_Matrix,scale = "none", border_color= "grey", number_color= "black",
                        fontsize_number=7,fontsize_row=10,fontsize_col=10,cellwidth=20,
                        cellheight=20,cluster_rows=F,cluster_cols=F,display_numbers = T,treeheight_row = 10,treeheight_col = 10,
                        color= mycol,show_rownames=T) 
ggsave(filename = "Normalized confusion matrix1.pdf", p1,device = 'pdf', width = 10, height = 10, units = 'cm')

# test
test <- recipe(Group~.,test_expr) %>%           
  step_dummy(all_nominal_predictors()) %>%           
  prep() %>%           
  bake(new_data=NULL) %>%          
  as.h2o() 
test$Group <- h2o.asfactor(test$Group)
y <- "Group" 
x <- setdiff(names(test),y)  
pred <- h2o.predict(best, test)
## 
Confusion_Matrix <- as.matrix(table(cbind(as.data.frame(pred$predict), as.data.frame(test$Group))))
Confusion_Matrix <- prop.table(Confusion_Matrix, margin = 1)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap::pheatmap(Confusion_Matrix,scale = "none", border_color= "grey", number_color= "black",
                        fontsize_number=7,fontsize_row=10,fontsize_col=10,cellwidth=20,
                        cellheight=20,cluster_rows=F,cluster_cols=F,display_numbers = T,treeheight_row = 10,treeheight_col = 10,
                        color= mycol,show_rownames=T) 
ggsave(filename = "Normalized confusion matrix2.pdf", p1,device = 'pdf', width = 10, height = 10, units = 'cm')


# ROC-AUC
h2o.auc(perf)
plot(perf,type="roc")
h2o.varimp_plot(best)
pdf("permutation_importance.pdf", width = 5,height = 5)
h2o.permutation_importance_plot(best,rec,num_of_features=5)
dev.off()

################################################################################
# Figure 8G scFEA-group analysis
################################################################################
library(ggpubr)
predME <- colnames(predFlux)
table(scRNA$Sex,scRNA$Organ)
table(scRNA$Sex)
current.cluster.ids <- c("F","female","M","male","NA")
new.cluster.ids <- c("F","F","M","M","NA") 
scRNA$Sex <- plyr::mapvalues(x = as.character(scRNA$Sex), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA$Sex)
Idents(scRNA) <- "Sex"
scRNA = subset(scRNA,idents = "NA", invert = TRUE)
scRNA$Sex <- as.factor(as.character(scRNA$Sex))
table(scRNA$Sex,scRNA$Organ)
scRNA = scRNA[,scRNA$Organ %in% c("Aorta","Bladder","Brain","Heart","Intestines","Kidney","Liver","Lung","Muscle","Pancreas","PIT","Retina","SAT","SC","Skin","VAT")]
scRNA$Organ <- as.factor(as.character(scRNA$Organ))
table(scRNA$Sex,scRNA$Organ)

{
{
#
data <- FetchData(scRNA, vars = c("Sex",predME[1],"Organ"))
p_ME <- compare_means(AMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME1 <- data.frame(Organ=p_ME$Organ,AMP=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Sex",predME[2],"Organ"))
p_ME <- compare_means(Pyruvate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME2 <- data.frame(Organ=p_ME$Organ,Pyruvate=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Sex",predME[3],"Organ"))
p_ME <- compare_means(Acetyl.CoA~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME3 <- data.frame(Organ=p_ME$Organ,Acetyl.CoA=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Sex",predME[4],"Organ"))
p_ME <- compare_means(Glutamate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME4 <- data.frame(Organ=p_ME$Organ,Glutamate=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Sex",predME[5],"Organ"))
p_ME <- compare_means(X2OG~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME5 <- data.frame(Organ=p_ME$Organ,X2OG=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Sex",predME[6],"Organ"))
p_ME <- compare_means(Oxaloacetate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME6 <- data.frame(Organ=p_ME$Organ,Oxaloacetate=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Sex",predME[7],"Organ"))
p_ME <- compare_means(Glycine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME7 <- data.frame(Organ=p_ME$Organ,Glycine=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Sex",predME[8],"Organ"))
p_ME <- compare_means(Succinate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME8 <- data.frame(Organ=p_ME$Organ,Succinate=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Sex",predME[9],"Organ"))
p_ME <- compare_means(UDP.N.acetylglucosamine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME9 <- data.frame(Organ=p_ME$Organ,UDP.N.acetylglucosamine=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Sex",predME[10],"Organ"))
p_ME <- compare_means(lysine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME10 <- data.frame(Organ=p_ME$Organ,lysine=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[11],"Organ"))
  p_ME <- compare_means(Aspartate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME11 <- data.frame(Organ=p_ME$Organ,Aspartate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[12],"Organ"))
  p_ME <- compare_means(Glutathione~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME12 <- data.frame(Organ=p_ME$Organ,Glutathione=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[13],"Organ"))
  p_ME <- compare_means(Arginine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME13 <- data.frame(Organ=p_ME$Organ,Arginine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[14],"Organ"))
  p_ME <- compare_means(Glutamine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME14 <- data.frame(Organ=p_ME$Organ,Glutamine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[15],"Organ"))
  p_ME <- compare_means(Serine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME15 <- data.frame(Organ=p_ME$Organ,Serine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[16],"Organ"))
  p_ME <- compare_means(Methionine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME16 <- data.frame(Organ=p_ME$Organ,Methionine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[17],"Organ"))
  p_ME <- compare_means(Ornithine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME17 <- data.frame(Organ=p_ME$Organ,Ornithine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[18],"Organ"))
  p_ME <- compare_means(Phenylalanine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME18 <- data.frame(Organ=p_ME$Organ,Phenylalanine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[19],"Organ"))
  p_ME <- compare_means(Tyrosine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME19 <- data.frame(Organ=p_ME$Organ,Tyrosine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[20],"Organ"))
  p_ME <- compare_means(Succinyl.CoA~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME20 <- data.frame(Organ=p_ME$Organ,Succinyl.CoA=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[21],"Organ"))
  p_ME <- compare_means(Cysteine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME21 <- data.frame(Organ=p_ME$Organ,Cysteine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[22],"Organ"))
  p_ME <- compare_means(B.Alanine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME22 <- data.frame(Organ=p_ME$Organ,B.Alanine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[23],"Organ"))
  p_ME <- compare_means(Propanoyl.CoA~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME23 <- data.frame(Organ=p_ME$Organ,Propanoyl.CoA=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[24],"Organ"))
  p_ME <- compare_means(Glucose.1.phosphate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME24 <- data.frame(Organ=p_ME$Organ,Glucose.1.phosphate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[25],"Organ"))
  p_ME <- compare_means(UMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME25 <- data.frame(Organ=p_ME$Organ,UMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[26],"Organ"))
  p_ME <- compare_means(Dolichyl.phosphate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME26 <- data.frame(Organ=p_ME$Organ,Dolichyl.phosphate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[27],"Organ"))
  p_ME <- compare_means(CDP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME27 <- data.frame(Organ=p_ME$Organ,CDP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[28],"Organ"))
  p_ME <- compare_means(Choline~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME28 <- data.frame(Organ=p_ME$Organ,Choline=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[29],"Organ"))
  p_ME <- compare_means(G3P~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME29 <- data.frame(Organ=p_ME$Organ,G3P=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[30],"Organ"))
  p_ME <- compare_means(PRPP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME30 <- data.frame(Organ=p_ME$Organ,PRPP=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[31],"Organ"))
  p_ME <- compare_means(Fumarate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME31 <- data.frame(Organ=p_ME$Organ,Fumarate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[32],"Organ"))
  p_ME <- compare_means(Leucine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME32 <- data.frame(Organ=p_ME$Organ,Leucine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[33],"Organ"))
  p_ME <- compare_means(IMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME33 <- data.frame(Organ=p_ME$Organ,IMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[34],"Organ"))
  p_ME <- compare_means(Putrescine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME34 <- data.frame(Organ=p_ME$Organ,Putrescine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[35],"Organ"))
  p_ME <- compare_means(Histidine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME35 <- data.frame(Organ=p_ME$Organ,Histidine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[36],"Organ"))
  p_ME <- compare_means(GMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME36 <- data.frame(Organ=p_ME$Organ,GMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[37],"Organ"))
  p_ME <- compare_means(Proline~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME37 <- data.frame(Organ=p_ME$Organ,Proline=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[38],"Organ"))
  p_ME <- compare_means(Malate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME38 <- data.frame(Organ=p_ME$Organ,Malate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[39],"Organ"))
  p_ME <- compare_means(Citrate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME39 <- data.frame(Organ=p_ME$Organ,Citrate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[40],"Organ"))
  p_ME <- compare_means(UDP.glucuronic.acid~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME40 <- data.frame(Organ=p_ME$Organ,UDP.glucuronic.acid=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[41],"Organ"))
  p_ME <- compare_means(Valine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME41 <- data.frame(Organ=p_ME$Organ,Valine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[42],"Organ"))
  p_ME <- compare_means(Cholesterol~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME42 <- data.frame(Organ=p_ME$Organ,Cholesterol=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[43],"Organ"))
  p_ME <- compare_means(Threonine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME43 <- data.frame(Organ=p_ME$Organ,Threonine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[44],"Organ"))
  p_ME <- compare_means(X3PD~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME44 <- data.frame(Organ=p_ME$Organ,X3PD=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[45],"Organ"))
  p_ME <- compare_means(dCMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME45 <- data.frame(Organ=p_ME$Organ,dCMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[46],"Organ"))
  p_ME <- compare_means(Lactate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME46 <- data.frame(Organ=p_ME$Organ,Lactate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[47],"Organ"))
  p_ME <- compare_means(Hypoxanthine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME47 <- data.frame(Organ=p_ME$Organ,Hypoxanthine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[48],"Organ"))
  p_ME <- compare_means(Glucose~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME48 <- data.frame(Organ=p_ME$Organ,Glucose=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[49],"Organ"))
  p_ME <- compare_means(Citrulline~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME49 <- data.frame(Organ=p_ME$Organ,Citrulline=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[50],"Organ"))
  p_ME <- compare_means(GABA~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME50 <- data.frame(Organ=p_ME$Organ,GABA=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[51],"Organ"))
  p_ME <- compare_means(dTMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME51 <- data.frame(Organ=p_ME$Organ,dTMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[52],"Organ"))
  p_ME <- compare_means(dUMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME52 <- data.frame(Organ=p_ME$Organ,dUMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[53],"Organ"))
  p_ME <- compare_means(Xanthine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME53 <- data.frame(Organ=p_ME$Organ,Xanthine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[54],"Organ"))
  p_ME <- compare_means(Chondroitin~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME54 <- data.frame(Organ=p_ME$Organ,Chondroitin=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[55],"Organ"))
  p_ME <- compare_means(Isoleucine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME55 <- data.frame(Organ=p_ME$Organ,Isoleucine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[56],"Organ"))
  p_ME <- compare_means(X.E.E..Farnesyl.PP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME56 <- data.frame(Organ=p_ME$Organ,X.E.E..Farnesyl.PP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[57],"Organ"))
  p_ME <- compare_means(Deoxyadenosine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME57 <- data.frame(Organ=p_ME$Organ,Deoxyadenosine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[58],"Organ"))
  p_ME <- compare_means(XMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME58 <- data.frame(Organ=p_ME$Organ,XMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[59],"Organ"))
  p_ME <- compare_means(G6P~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME59 <- data.frame(Organ=p_ME$Organ,G6P=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[60],"Organ"))
  p_ME <- compare_means(dCDP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME60 <- data.frame(Organ=p_ME$Organ,dCDP=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[61],"Organ"))
  p_ME <- compare_means(Spermine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME61 <- data.frame(Organ=p_ME$Organ,Spermine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[62],"Organ"))
  p_ME <- compare_means(Argininosuccinate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME62 <- data.frame(Organ=p_ME$Organ,Argininosuccinate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[63],"Organ"))
  p_ME <- compare_means(Dolichyl.phosphate.D.mannose~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME63 <- data.frame(Organ=p_ME$Organ,Dolichyl.phosphate.D.mannose=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[64],"Organ"))
  p_ME <- compare_means(AICAR~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME64 <- data.frame(Organ=p_ME$Organ,AICAR=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[65],"Organ"))
  p_ME <- compare_means(Fatty.Acid~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME65 <- data.frame(Organ=p_ME$Organ,Fatty.Acid=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[66],"Organ"))
  p_ME <- compare_means(Pyrimidine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME66 <- data.frame(Organ=p_ME$Organ,Pyrimidine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[67],"Organ"))
  p_ME <- compare_means(X.Glc.3..GlcNAc.2..Man.9..Asn.1~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME67 <- data.frame(Organ=p_ME$Organ,X.Glc.3..GlcNAc.2..Man.9..Asn.1=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[68],"Organ"))
  p_ME <- compare_means(X.GlcNAc.4..Man.3..Asn.1~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME68 <- data.frame(Organ=p_ME$Organ,X.GlcNAc.4..Man.3..Asn.1=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[69],"Organ"))
  p_ME <- compare_means(X.Gal.2..GlcA.1..Xyl.1..Ser.1~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME69 <- data.frame(Organ=p_ME$Organ,X.Gal.2..GlcA.1..Xyl.1..Ser.1=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Sex",predME[70],"Organ"))
  p_ME <- compare_means(X.Gal.1..GlcNAc.1..Man.1..Ser.Thr.1~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME70 <- data.frame(Organ=p_ME$Organ,X.Gal.1..GlcNAc.1..Man.1..Ser.Thr.1=p_ME$p.adj)
}

{
p_ME_all <- merge(p_ME1,  p_ME2, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME3, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME4, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME5, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME6, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME7, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME8, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME9, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME10, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME11, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME12, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME13, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME14, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME15, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME16, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME17, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME18, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME19, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME20, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME21, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME22, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME23, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME24, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME25, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME26, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME27, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME28, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME29, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME30, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME31, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME32, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME33, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME34, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME35, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME36, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME37, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME38, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME39, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME40, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME41, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME42, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME43, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME44, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME45, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME46, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME47, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME48, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME49, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME50, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME51, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME52, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME53, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME54, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME55, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME56, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME57, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME58, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME59, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME60, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME61, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME62, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME63, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME64, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME65, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME66, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME67, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME68, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME69, by="Organ",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME70, by="Organ",all=TRUE)
}

head(p_ME_all)

write.csv(p_ME_all,"p_ME_all.csv", na="0",row.names = F) 
}

{
  {
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[1],"Organ"))
    p_ME <- compare_means(AMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME1 <- data.frame(Organ=p_ME$Organ,AMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[2],"Organ"))
    p_ME <- compare_means(Pyruvate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME2 <- data.frame(Organ=p_ME$Organ,Pyruvate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[3],"Organ"))
    p_ME <- compare_means(Acetyl.CoA~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME3 <- data.frame(Organ=p_ME$Organ,Acetyl.CoA=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[4],"Organ"))
    p_ME <- compare_means(Glutamate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME4 <- data.frame(Organ=p_ME$Organ,Glutamate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[5],"Organ"))
    p_ME <- compare_means(X2OG~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME5 <- data.frame(Organ=p_ME$Organ,X2OG=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[6],"Organ"))
    p_ME <- compare_means(Oxaloacetate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME6 <- data.frame(Organ=p_ME$Organ,Oxaloacetate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[7],"Organ"))
    p_ME <- compare_means(Glycine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME7 <- data.frame(Organ=p_ME$Organ,Glycine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[8],"Organ"))
    p_ME <- compare_means(Succinate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME8 <- data.frame(Organ=p_ME$Organ,Succinate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[9],"Organ"))
    p_ME <- compare_means(UDP.N.acetylglucosamine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME9 <- data.frame(Organ=p_ME$Organ,UDP.N.acetylglucosamine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[10],"Organ"))
    p_ME <- compare_means(lysine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME10 <- data.frame(Organ=p_ME$Organ,lysine=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[11],"Organ"))
    p_ME <- compare_means(Aspartate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME11 <- data.frame(Organ=p_ME$Organ,Aspartate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[12],"Organ"))
    p_ME <- compare_means(Glutathione~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME12 <- data.frame(Organ=p_ME$Organ,Glutathione=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[13],"Organ"))
    p_ME <- compare_means(Arginine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME13 <- data.frame(Organ=p_ME$Organ,Arginine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[14],"Organ"))
    p_ME <- compare_means(Glutamine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME14 <- data.frame(Organ=p_ME$Organ,Glutamine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[15],"Organ"))
    p_ME <- compare_means(Serine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME15 <- data.frame(Organ=p_ME$Organ,Serine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[16],"Organ"))
    p_ME <- compare_means(Methionine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME16 <- data.frame(Organ=p_ME$Organ,Methionine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[17],"Organ"))
    p_ME <- compare_means(Ornithine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME17 <- data.frame(Organ=p_ME$Organ,Ornithine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[18],"Organ"))
    p_ME <- compare_means(Phenylalanine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME18 <- data.frame(Organ=p_ME$Organ,Phenylalanine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[19],"Organ"))
    p_ME <- compare_means(Tyrosine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME19 <- data.frame(Organ=p_ME$Organ,Tyrosine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[20],"Organ"))
    p_ME <- compare_means(Succinyl.CoA~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME20 <- data.frame(Organ=p_ME$Organ,Succinyl.CoA=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[21],"Organ"))
    p_ME <- compare_means(Cysteine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME21 <- data.frame(Organ=p_ME$Organ,Cysteine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[22],"Organ"))
    p_ME <- compare_means(B.Alanine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME22 <- data.frame(Organ=p_ME$Organ,B.Alanine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[23],"Organ"))
    p_ME <- compare_means(Propanoyl.CoA~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME23 <- data.frame(Organ=p_ME$Organ,Propanoyl.CoA=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[24],"Organ"))
    p_ME <- compare_means(Glucose.1.phosphate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME24 <- data.frame(Organ=p_ME$Organ,Glucose.1.phosphate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[25],"Organ"))
    p_ME <- compare_means(UMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME25 <- data.frame(Organ=p_ME$Organ,UMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[26],"Organ"))
    p_ME <- compare_means(Dolichyl.phosphate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME26 <- data.frame(Organ=p_ME$Organ,Dolichyl.phosphate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[27],"Organ"))
    p_ME <- compare_means(CDP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME27 <- data.frame(Organ=p_ME$Organ,CDP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[28],"Organ"))
    p_ME <- compare_means(Choline~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME28 <- data.frame(Organ=p_ME$Organ,Choline=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[29],"Organ"))
    p_ME <- compare_means(G3P~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME29 <- data.frame(Organ=p_ME$Organ,G3P=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[30],"Organ"))
    p_ME <- compare_means(PRPP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME30 <- data.frame(Organ=p_ME$Organ,PRPP=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[31],"Organ"))
    p_ME <- compare_means(Fumarate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME31 <- data.frame(Organ=p_ME$Organ,Fumarate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[32],"Organ"))
    p_ME <- compare_means(Leucine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME32 <- data.frame(Organ=p_ME$Organ,Leucine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[33],"Organ"))
    p_ME <- compare_means(IMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME33 <- data.frame(Organ=p_ME$Organ,IMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[34],"Organ"))
    p_ME <- compare_means(Putrescine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME34 <- data.frame(Organ=p_ME$Organ,Putrescine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[35],"Organ"))
    p_ME <- compare_means(Histidine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME35 <- data.frame(Organ=p_ME$Organ,Histidine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[36],"Organ"))
    p_ME <- compare_means(GMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME36 <- data.frame(Organ=p_ME$Organ,GMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[37],"Organ"))
    p_ME <- compare_means(Proline~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME37 <- data.frame(Organ=p_ME$Organ,Proline=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[38],"Organ"))
    p_ME <- compare_means(Malate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME38 <- data.frame(Organ=p_ME$Organ,Malate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[39],"Organ"))
    p_ME <- compare_means(Citrate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME39 <- data.frame(Organ=p_ME$Organ,Citrate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[40],"Organ"))
    p_ME <- compare_means(UDP.glucuronic.acid~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME40 <- data.frame(Organ=p_ME$Organ,UDP.glucuronic.acid=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[41],"Organ"))
    p_ME <- compare_means(Valine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME41 <- data.frame(Organ=p_ME$Organ,Valine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[42],"Organ"))
    p_ME <- compare_means(Cholesterol~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME42 <- data.frame(Organ=p_ME$Organ,Cholesterol=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[43],"Organ"))
    p_ME <- compare_means(Threonine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME43 <- data.frame(Organ=p_ME$Organ,Threonine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[44],"Organ"))
    p_ME <- compare_means(X3PD~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME44 <- data.frame(Organ=p_ME$Organ,X3PD=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[45],"Organ"))
    p_ME <- compare_means(dCMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME45 <- data.frame(Organ=p_ME$Organ,dCMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[46],"Organ"))
    p_ME <- compare_means(Lactate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME46 <- data.frame(Organ=p_ME$Organ,Lactate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[47],"Organ"))
    p_ME <- compare_means(Hypoxanthine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME47 <- data.frame(Organ=p_ME$Organ,Hypoxanthine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[48],"Organ"))
    p_ME <- compare_means(Glucose~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME48 <- data.frame(Organ=p_ME$Organ,Glucose=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[49],"Organ"))
    p_ME <- compare_means(Citrulline~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME49 <- data.frame(Organ=p_ME$Organ,Citrulline=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[50],"Organ"))
    p_ME <- compare_means(GABA~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME50 <- data.frame(Organ=p_ME$Organ,GABA=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[51],"Organ"))
    p_ME <- compare_means(dTMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME51 <- data.frame(Organ=p_ME$Organ,dTMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[52],"Organ"))
    p_ME <- compare_means(dUMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME52 <- data.frame(Organ=p_ME$Organ,dUMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[53],"Organ"))
    p_ME <- compare_means(Xanthine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME53 <- data.frame(Organ=p_ME$Organ,Xanthine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[54],"Organ"))
    p_ME <- compare_means(Chondroitin~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME54 <- data.frame(Organ=p_ME$Organ,Chondroitin=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[55],"Organ"))
    p_ME <- compare_means(Isoleucine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME55 <- data.frame(Organ=p_ME$Organ,Isoleucine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[56],"Organ"))
    p_ME <- compare_means(X.E.E..Farnesyl.PP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME56 <- data.frame(Organ=p_ME$Organ,X.E.E..Farnesyl.PP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[57],"Organ"))
    p_ME <- compare_means(Deoxyadenosine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME57 <- data.frame(Organ=p_ME$Organ,Deoxyadenosine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[58],"Organ"))
    p_ME <- compare_means(XMP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME58 <- data.frame(Organ=p_ME$Organ,XMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[59],"Organ"))
    p_ME <- compare_means(G6P~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME59 <- data.frame(Organ=p_ME$Organ,G6P=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[60],"Organ"))
    p_ME <- compare_means(dCDP~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME60 <- data.frame(Organ=p_ME$Organ,dCDP=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[61],"Organ"))
    p_ME <- compare_means(Spermine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME61 <- data.frame(Organ=p_ME$Organ,Spermine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[62],"Organ"))
    p_ME <- compare_means(Argininosuccinate~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME62 <- data.frame(Organ=p_ME$Organ,Argininosuccinate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[63],"Organ"))
    p_ME <- compare_means(Dolichyl.phosphate.D.mannose~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME63 <- data.frame(Organ=p_ME$Organ,Dolichyl.phosphate.D.mannose=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[64],"Organ"))
    p_ME <- compare_means(AICAR~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME64 <- data.frame(Organ=p_ME$Organ,AICAR=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[65],"Organ"))
    p_ME <- compare_means(Fatty.Acid~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME65 <- data.frame(Organ=p_ME$Organ,Fatty.Acid=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[66],"Organ"))
    p_ME <- compare_means(Pyrimidine~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME66 <- data.frame(Organ=p_ME$Organ,Pyrimidine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[67],"Organ"))
    p_ME <- compare_means(X.Glc.3..GlcNAc.2..Man.9..Asn.1~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME67 <- data.frame(Organ=p_ME$Organ,X.Glc.3..GlcNAc.2..Man.9..Asn.1=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[68],"Organ"))
    p_ME <- compare_means(X.GlcNAc.4..Man.3..Asn.1~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME68 <- data.frame(Organ=p_ME$Organ,X.GlcNAc.4..Man.3..Asn.1=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[69],"Organ"))
    p_ME <- compare_means(X.Gal.2..GlcA.1..Xyl.1..Ser.1~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME69 <- data.frame(Organ=p_ME$Organ,X.Gal.2..GlcA.1..Xyl.1..Ser.1=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Sex",predME[70],"Organ"))
    p_ME <- compare_means(X.Gal.1..GlcNAc.1..Man.1..Ser.Thr.1~Sex, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME70 <- data.frame(Organ=p_ME$Organ,X.Gal.1..GlcNAc.1..Man.1..Ser.Thr.1=-log10(p_ME$p.adj))
  }
  
  {
    p_ME_all <- merge(p_ME1,  p_ME2, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME3, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME4, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME5, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME6, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME7, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME8, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME9, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME10, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME11, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME12, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME13, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME14, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME15, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME16, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME17, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME18, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME19, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME20, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME21, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME22, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME23, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME24, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME25, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME26, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME27, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME28, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME29, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME30, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME31, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME32, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME33, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME34, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME35, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME36, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME37, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME38, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME39, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME40, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME41, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME42, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME43, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME44, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME45, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME46, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME47, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME48, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME49, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME50, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME51, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME52, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME53, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME54, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME55, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME56, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME57, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME58, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME59, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME60, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME61, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME62, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME63, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME64, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME65, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME66, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME67, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME68, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME69, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME70, by="Organ",all=TRUE)
  }
  
  head(p_ME_all)
  
  write.csv(p_ME_all,"logp_ME_all.csv", na="0",row.names = F)
}

r2 <-as.matrix(read.csv('logp_ME_all.csv', header = T, row.names = 1))
r2[is.infinite(r2)] <- 300
r2[r2 >= 10] <- 10
r3 <- read.csv('scFEA_all_diff.csv', header = T, row.names = 1)
r3 <- r3*-1
r3[r3 > 0] <- 1
r3[r3 < 0] <- -1
r4 <- r3*r2
p2 <- read.csv('p_ME_all.csv', header = T, row.names = 1)
p2[p2 < 0.001] <- "***"
p2[p2 < 0.01 & p2 >= 0.001] <- "**"
p2[p2 < 0.05 & p2 >= 0.01] <- "*"
p2[p2 >= 0.05] <- ""
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap::pheatmap(r4,scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=10,fontsize_row=8,fontsize_col=9,cellwidth=14,
              cellheight=12,cluster_rows=T,cluster_cols=T,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "p_ME_all1.pdf", p1,device = 'pdf', width = 50, height = 20, units = 'cm')



