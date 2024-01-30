 
################################################################################
# Figure 7B DEGs
################################################################################
# MAST was used to calculate senescence-related differential genes in ECs from different tissues
library(MAST);library(Seurat);library(dplyr)
load("Brain_Endo1.RData")
Idents(scRNA)="celltype"
table(Idents(scRNA))
table(scRNA$Group)
scRNA1 = scRNA
fData = data.frame(symbolid=rownames(scRNA1),primerid=rownames(scRNA1))
rownames(fData)=fData$symbolid
cData = scRNA1@meta.data
cData$wellKey <- rownames(cData)
sca = FromMatrix(as.matrix(scRNA1@assays$RNA@data), cData = cData,fData = fData)
rm(scRNA1)
gc()
dim(sca)
table(colData(sca)$Group)
cond<-factor(colData(sca)$Group)
cond<-relevel(cond,"Adult")
colData(sca)$condition<-cond

colData(sca)$Sex=factor(colData(sca)$Sex)
colData(sca)$Age=as.numeric(colData(sca)$Age)
colData(sca)$percent.mt=as.numeric(colData(sca)$percent.mt)
colData(sca)$percent.rb=as.numeric(colData(sca)$percent.rb)
colData(sca)$nCount_RNA=as.numeric(colData(sca)$nCount_RNA)
#colData(sca)$datasets=factor(colData(sca)$datasets)

##
zlmCond <- zlm(~condition + nCount_RNA + percent.mt + percent.rb + Sex, sca, method="bayesglm", ebayes=TRUE)
summaryCond <- summary(zlmCond,doLRT='conditionAging')
summaryDt <- summaryCond$datatable
levels(summaryDt$contrast)

#
df_pval = summaryDt %>% 
  dplyr::filter(contrast=='conditionAging') %>% 
  dplyr::filter(component=='H') %>% 
  dplyr::select(primerid, `Pr(>Chisq)`)

df_logfc = summaryDt %>% 
  dplyr::filter(contrast=='conditionAging') %>% 
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
write.csv(df_stat,"Brain_Aging_Adult_DEGs.csv")

################################################################################
# Figure 7C Overlap of genes associated with aging from different tissues
################################################################################
library(dplyr);library(GeneOverlap);data(GeneOverlap)
#Select variables
variable_list=list("Group")
for (k in 1:length(variable_list)) {
  variable=variable_list[[k]]
  
  setwd("/Users/niuruize/Downloads/scRNA/EC/EC/7_aging")
  
  celltype = c("Aorta","Bladder","Brain","Heart","Intestinal","Kidney","Liver","Lung","Muscle","Pancreas","PIT","Retina","SAT","SC","Skin","Testis","Thyroid","Uterus","VAT")
  
  #List significantly upregulated genes
  
  Marker_genes_list<-list()
  for (i in 1:length(celltype)) {
    ct1=celltype[i]
    ct1=gsub("/","_",ct1)
    input_file=paste0(ct1,"_Aging_Adult_DEGs",".csv")
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
      
      input_file1=paste0(ct1_name,"_Aging_Adult_DEGs",".csv")
      input_file2=paste0(ct2_name,"_Aging_Adult_DEGs",".csv")
      
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
  
  setwd("/Users/niuruize/Downloads/scRNA/EC/EC/7_aging")
  output_file=paste0("EC_aging_up_P",".csv")
  write.csv(pvalue_adjusted_matrix_minusLog10,file=output_file)
  
  
  GeneOverlap_Matrix=getMatrix(gom.obj, name="odds.ratio")
  GeneOverlap_Matrix[is.infinite(GeneOverlap_Matrix)] <- 20
  odds_ratio_Matrix = as.data.frame(GeneOverlap_Matrix)
  output_file=paste0("EC_aging_up_OR",".csv")
  write.csv(odds_ratio_Matrix,file=output_file)
}

library(pheatmap)
library(ggplot2)
r2 <- read.csv('EC_aging_up_OR.csv',header=TRUE,row.names=1)
p2 <- read.csv('EC_aging_up_P.csv',header=TRUE,row.names=1)
r2[r2 >= 10] <- 10
p2[p2 <= 30] <- ""
p2[p2 >= 100] <- "**"
p2[p2 > 30] <- "*"
#mycol<-colorRampPalette(c( "#0f86a9", "white", "#ed8b10"))(200)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
#
p1 = pheatmap(r2,scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=14,fontsize_row=8,fontsize_col=9,cellwidth=15,
              cellheight=15,cluster_rows=T,cluster_cols=T,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "EC_aging_overlap_up1.pdf", p1,device = 'pdf', width = 15, height = 15, units = 'cm')
p1 = pheatmap(r2,scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=14,fontsize_row=8,fontsize_col=9,cellwidth=15,
              cellheight=15,cluster_rows=F,cluster_cols=F,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "EC_aging_overlap_up2.pdf", p1,device = 'pdf', width = 15, height = 15, units = 'cm')

################################################################################
# Figure 7E GO analysis
################################################################################
library("clusterProfiler");library("org.Hs.eg.db");library("enrichplot");library("ggplot2");library("org.Hs.eg.db")
##
#--------------------------------------------------------------------Dowm---------------------------------------------------
celltype = c("Aorta","Bladder","Brain","Heart","Intestinal","Kidney","Liver","Lung","Muscle","Pancreas","Retina","SAT","SC","Skin","Testis","Thyroid","Uterus","VAT")
for (i in 1:length(celltype)) {
  ct1=celltype[i]
  ct1=gsub("/","_",ct1)
  input_file=paste0("/Users/niuruize/Downloads/scRNA/EC/EC/7_aging/",ct1,"_Aging_Adult_DEGs",".csv")
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
  input_file=paste0("/Users/niuruize/Downloads/scRNA/EC/EC/7_aging/",ct1,"_Aging_Adult_DEGs",".csv")
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
# Figure 7F Virtual KO
################################################################################
# Extract expression matrix
library(scTenifoldKnk);library(Seurat)
load(Skin_Harmony.RData)
table(scRNA$Group)
scRNA_seq <- scRNA[,scRNA$Group %in% c("Aging")]
scRNA_seq$Group <- as.factor(as.character(scRNA_seq$Group))
X <- data.frame(scRNA_seq@assays$RNA@counts)
save(X ,file="Skin_Aging_KO_data.RData")
# Virtual KO
library(scTenifoldKnk);library(Seurat);library(Matrix)
args = commandArgs(trailingOnly=TRUE)
scQC <- function(X, mtThreshold = 0.1, minLSize = 1000){
  require(Matrix)
  if(class(X) == 'Seurat'){
    countMatrix <- X@assays$RNA@counts
  } else {
    countMatrix <- X
  }
  librarySize <- colSums(countMatrix)
  countMatrix <- countMatrix[,librarySize >= minLSize]
  librarySize <- colSums(countMatrix)
  mtGenes <- grep('^MT-',toupper(rownames(countMatrix)))
  nGenes <- colSums(countMatrix != 0)
  
  genesLM <- lm(nGenes~librarySize)
  genesLM <- as.data.frame(predict(genesLM, data.frame(librarySize), interval = 'prediction'))
  
  if(isTRUE(length(mtGenes) > 0)){
    mtCounts <- colSums(countMatrix[grep('^MT-',toupper(rownames(countMatrix))),])
    mtProportion <- mtCounts/librarySize
    mtLM <- lm(mtCounts~librarySize)
    mtLM <- as.data.frame(predict(mtLM, data.frame(librarySize), interval = 'prediction'))
    selectedCells <- ((mtCounts > mtLM$lwr) & (mtCounts < mtLM$upr) & (nGenes > genesLM$lwr) & (nGenes < genesLM$upr) & (mtProportion <= mtThreshold) & (librarySize < 2 * mean(librarySize)))
  } else {
    selectedCells <- ((nGenes > genesLM$lwr) & (nGenes < genesLM$upr) & (librarySize < 2 * mean(librarySize)))
  }
  selectedCells <- colnames(countMatrix)[selectedCells]
  if(class(X) == 'Seurat'){
    X <- subset(X, cells = selectedCells)
  } else {
    X <- countMatrix[,selectedCells]
  }
  return(X)
}

# Running scTenifoldKnk
load("~/aging/EC/KO/Skin_KO_counts.RData")
X <- scQC(X)
tempX <- X[,sample(seq_len(ncol(X)),146)]
tempX <- tempX[rowMeans(X != 0) > 0.05,]
oX <- scTenifoldKnk(countMatrix = tempX, gKO = 'JAM2',qc = TRUE, nCores = 16)
save(oX ,file="Skin_JAM2_counts.RData")
write.csv(oX$diffRegulation, file = paste0('JAM2_counts_Stability_',args[1], '.csv'), row.names = FALSE)
# visualization
library(fgsea)；library(ggplot2)；library(enrichR)；library(igraph)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotKO.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa2mmu_SYMBOL.R')

load('Skin_JAM2.RData')
write.csv(Skin_JAM2$diffRegulation, '../Results/drTREM2.csv')
dGenes <- Skin_JAM2$diffRegulation$gene[Skin_JAM2$diffRegulation$p.adj < 0.05]

#png('dr2_Skin_JAM2.png', width = 2000, height = 2000, res = 300)
pdf("dr2_Skin_JAM2.pdf", width = 5,height = 5)
plotDR(Skin_JAM2)
dev.off()

#png('ego2_Skin_JAM2.png', width = 3000, height = 3000, res = 300, bg = NA)
pdf("ego2_Skin_JAM2.pdf", width = 10,height = 10)
X <- Skin_JAM2
gKO <- 'JAM2'
q <- 0.995
gList <- unique(c(gKO, X$diffRegulation$gene[X$diffRegulation$p.adj < 0.05]))
sCluster <- as.matrix(X$tensorNetworks$WT[gList,gList])
koInfo <- sCluster[gKO,]
sCluster[abs(sCluster) <= quantile(abs(sCluster), q)] <- 0
sCluster[gKO,] <- koInfo
diag(sCluster) <- 0
sCluster <-  reshape2::melt(as.matrix(sCluster))
colnames(sCluster) <- c('from', 'to', 'W')
sCluster <- sCluster[sCluster$W != 0,]
netPlot <- graph_from_data_frame(sCluster, directed = TRUE)
dPlot <- centr_degree(netPlot)$res
W <- rep(1,nrow(sCluster))
sG   <- (names(V(netPlot))[dPlot > 1])[-1]
W[sCluster$from %in% sG] <- 0.2
W[sCluster$to %in% sG] <- 0.2
W[sCluster$from %in% gKO] <- 1
W[sCluster$from %in% gKO & sCluster$to %in% sG] <- 0.8
set.seed(1)
layPlot <- layout_with_fr(netPlot, weights = W)
dPlot <- (dPlot/max(dPlot))*20
dbs <- listEnrichrDbs()
E <- enrichr(gList, c("GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Biological_Process_2023","KEGG_2021_Human"))
E <- do.call(rbind.data.frame, E)
#E <- E[E$Adjusted.P.value < 0.05,]
E <- E[order(E$Adjusted.P.value),]
E$Term <- unlist(lapply(strsplit(E$Term,''), function(X){
  X[1] <- toupper(X[1])
  X <- paste0(X,collapse = '')
  X <- gsub('\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+','',X)
  X <- gsub("'s",'',X)
  X <- unlist(strsplit(X,','))[1]
  X <- gsub('[[:blank:]]$','',X)
  return(X)
}))
write.csv(E,"JAM2_EGO.csv") #保存结果
E <- E[E$Term %in% c('Antigen processing and presentation','MHC Class II Protein Complex Assembly','Positive Regulation Of Leukocyte Cell-Cell Adhesion','Lymphatic Endothelial Cell Differentiation','Regulation Of Endothelial Cell Migration','Regulation Of Endothelial Cell Proliferation','Regulation Of Vasculogenesis'),]
E <- E[c(1,2,4,6,8,10,12),]
tPlot <- strsplit(E$Genes, ';')
pPlot <- matrix(0,nrow = length(V(netPlot)), ncol = nrow(E))
rownames(pPlot) <- toupper(names(V(netPlot)))
for(i in seq_along(tPlot)){
  pPlot[unlist(tPlot[i]),i] <- 1
}
pPlot <- lapply(seq_len(nrow(pPlot)), function(X){as.vector(pPlot[X,])})
names(pPlot) <- names(V(netPlot))
tPlot <- unique(unlist(tPlot))
eGenes <- toupper(names(V(netPlot))) %in% tPlot
vColor <- rgb(195/255, 199/255, 198/255 ,0.3)
pieColors <- list(hcl.colors(nrow(E), palette = 'Zissou 1', alpha = 0.7))
par(mar=c(4,0,0,0), xpd = TRUE)
suppressWarnings(plot(netPlot,
                      layout = layPlot,
                      edge.arrow.size=.2,
                      vertex.label.color="black",
                      vertex.shape = ifelse(eGenes,'pie','circle'),
                      vertex.pie = pPlot,
                      vertex.size = 10+dPlot,
                      vertex.pie.color=pieColors,
                      vertex.label.family="Times",
                      vertex.label.font=ifelse(eGenes,2,1),
                      edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue'),
                      edge.curved = ifelse(W == 0.2, 0, 0.1),
                      vertex.color = vColor,
                      vertex.frame.color = NA))
sigLevel <- formatC(E$Adjusted.P.value, digits = 2, format = 'g', width = 0, drop0trailing = TRUE)
gSetNames <- lengths(strsplit(E$Genes, ';'))
gSetNames <- paste0('(', gSetNames,') ', E$Term, ' FDR = ', sigLevel)
legend(x = -1.05, y = -1.05, legend = gSetNames, bty = 'n', ncol = 2, cex = 1, col = unlist(pieColors), pch = 16)
dev.off()


zSkin_JAM2 <- Skin_JAM2$diffRegulation$Z
names(zSkin_JAM2) <- toupper(Skin_JAM2$diffRegulation$gene)

MGI <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human')
set.seed(1)
E <- fgseaMultilevel(MGI, zSkin_JAM2)

#png('gsea1_Skin_JAM2.png', width = 1000, height = 1000, res = 300)
pdf("gsea1_Skin_JAM2.pdf", width = 5,height = 5)
gSet <- 'Cell adhesion molecules (CAMs)'
plotEnrichment(MGI[[gSet]], zSkin_JAM2) +
  labs(
    title = 'Cell adhesion molecules (CAMs)',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea2_Skin_JAM2.png', width = 1000, height = 1000, res = 300)
pdf("gsea2_Skin_JAM2.pdf", width = 5,height = 5)
gSet <- 'Antigen processing and presentation'
plotEnrichment(MGI[[gSet]], zSkin_JAM2) +
  labs(
    title = 'Antigen processing and presentation',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea3_Skin_JAM2.png', width = 1000, height = 1000, res = 300)
pdf("gsea3_Skin_JAM2.pdf", width = 5,height = 5)
gSet <- 'Leukocyte transendothelial migration'
plotEnrichment(MGI[[gSet]], zSkin_JAM2) +
  labs(
    title = 'Leukocyte transendothelial migration',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score')  + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

################################################################################
# Figure 7G scFEA-group analysis
################################################################################
library(ggpubr)
predME <- colnames(predFlux)
table(scRNA$Group,scRNA$Organ)
table(scRNA$Group)
scRNA = scRNA[,scRNA$Organ %in% c("Aorta","Bladder","Brain","Heart","Intestines","Kidney","Liver","Lung","Muscle","Pancreas","PIT","Retina","SAT","SC","Skin","Testis","Thyroid","Uterus","VAT")]
scRNA$Organ <- as.factor(as.character(scRNA$Organ))
table(scRNA$Group,scRNA$Organ)
scRNA$Organ_Age <- paste(scRNA$Organ, scRNA$Age, sep = "_")
asd=data.frame(table(scRNA$Organ_Age))
write.csv(asd,"Organ_Age.csv", row.names = F) 
cluster_celltype <-  data.frame(readxl::read_xlsx("Organ_Age.xlsx"))
current.cluster.ids <- cluster_celltype$Var1
new.cluster.ids <- cluster_celltype$Group
scRNA$Group <- plyr::mapvalues(x = as.character(scRNA$Organ_Age), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA$Group,scRNA$Organ)
scRNA = scRNA[,scRNA$Group %in% c("Adult","Aging")]
scRNA$Group <- as.factor(as.character(scRNA$Group))
table(scRNA$Group,scRNA$Organ)

{
{
#
data <- FetchData(scRNA, vars = c("Group",predME[1],"Organ"))
p_ME <- compare_means(AMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME1 <- data.frame(Organ=p_ME$Organ,AMP=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[2],"Organ"))
p_ME <- compare_means(Pyruvate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME2 <- data.frame(Organ=p_ME$Organ,Pyruvate=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[3],"Organ"))
p_ME <- compare_means(Acetyl.CoA~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME3 <- data.frame(Organ=p_ME$Organ,Acetyl.CoA=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[4],"Organ"))
p_ME <- compare_means(Glutamate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME4 <- data.frame(Organ=p_ME$Organ,Glutamate=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[5],"Organ"))
p_ME <- compare_means(X2OG~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME5 <- data.frame(Organ=p_ME$Organ,X2OG=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[6],"Organ"))
p_ME <- compare_means(Oxaloacetate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME6 <- data.frame(Organ=p_ME$Organ,Oxaloacetate=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[7],"Organ"))
p_ME <- compare_means(Glycine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME7 <- data.frame(Organ=p_ME$Organ,Glycine=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[8],"Organ"))
p_ME <- compare_means(Succinate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME8 <- data.frame(Organ=p_ME$Organ,Succinate=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[9],"Organ"))
p_ME <- compare_means(UDP.N.acetylglucosamine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME9 <- data.frame(Organ=p_ME$Organ,UDP.N.acetylglucosamine=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[10],"Organ"))
p_ME <- compare_means(lysine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
p_ME10 <- data.frame(Organ=p_ME$Organ,lysine=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Group",predME[11],"Organ"))
  p_ME <- compare_means(Aspartate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME11 <- data.frame(Organ=p_ME$Organ,Aspartate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[12],"Organ"))
  p_ME <- compare_means(Glutathione~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME12 <- data.frame(Organ=p_ME$Organ,Glutathione=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[13],"Organ"))
  p_ME <- compare_means(Arginine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME13 <- data.frame(Organ=p_ME$Organ,Arginine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[14],"Organ"))
  p_ME <- compare_means(Glutamine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME14 <- data.frame(Organ=p_ME$Organ,Glutamine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[15],"Organ"))
  p_ME <- compare_means(Serine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME15 <- data.frame(Organ=p_ME$Organ,Serine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[16],"Organ"))
  p_ME <- compare_means(Methionine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME16 <- data.frame(Organ=p_ME$Organ,Methionine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[17],"Organ"))
  p_ME <- compare_means(Ornithine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME17 <- data.frame(Organ=p_ME$Organ,Ornithine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[18],"Organ"))
  p_ME <- compare_means(Phenylalanine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME18 <- data.frame(Organ=p_ME$Organ,Phenylalanine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[19],"Organ"))
  p_ME <- compare_means(Tyrosine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME19 <- data.frame(Organ=p_ME$Organ,Tyrosine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[20],"Organ"))
  p_ME <- compare_means(Succinyl.CoA~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME20 <- data.frame(Organ=p_ME$Organ,Succinyl.CoA=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Group",predME[21],"Organ"))
  p_ME <- compare_means(Cysteine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME21 <- data.frame(Organ=p_ME$Organ,Cysteine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[22],"Organ"))
  p_ME <- compare_means(B.Alanine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME22 <- data.frame(Organ=p_ME$Organ,B.Alanine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[23],"Organ"))
  p_ME <- compare_means(Propanoyl.CoA~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME23 <- data.frame(Organ=p_ME$Organ,Propanoyl.CoA=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[24],"Organ"))
  p_ME <- compare_means(Glucose.1.phosphate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME24 <- data.frame(Organ=p_ME$Organ,Glucose.1.phosphate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[25],"Organ"))
  p_ME <- compare_means(UMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME25 <- data.frame(Organ=p_ME$Organ,UMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[26],"Organ"))
  p_ME <- compare_means(Dolichyl.phosphate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME26 <- data.frame(Organ=p_ME$Organ,Dolichyl.phosphate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[27],"Organ"))
  p_ME <- compare_means(CDP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME27 <- data.frame(Organ=p_ME$Organ,CDP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[28],"Organ"))
  p_ME <- compare_means(Choline~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME28 <- data.frame(Organ=p_ME$Organ,Choline=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[29],"Organ"))
  p_ME <- compare_means(G3P~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME29 <- data.frame(Organ=p_ME$Organ,G3P=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[30],"Organ"))
  p_ME <- compare_means(PRPP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME30 <- data.frame(Organ=p_ME$Organ,PRPP=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Group",predME[31],"Organ"))
  p_ME <- compare_means(Fumarate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME31 <- data.frame(Organ=p_ME$Organ,Fumarate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[32],"Organ"))
  p_ME <- compare_means(Leucine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME32 <- data.frame(Organ=p_ME$Organ,Leucine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[33],"Organ"))
  p_ME <- compare_means(IMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME33 <- data.frame(Organ=p_ME$Organ,IMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[34],"Organ"))
  p_ME <- compare_means(Putrescine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME34 <- data.frame(Organ=p_ME$Organ,Putrescine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[35],"Organ"))
  p_ME <- compare_means(Histidine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME35 <- data.frame(Organ=p_ME$Organ,Histidine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[36],"Organ"))
  p_ME <- compare_means(GMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME36 <- data.frame(Organ=p_ME$Organ,GMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[37],"Organ"))
  p_ME <- compare_means(Proline~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME37 <- data.frame(Organ=p_ME$Organ,Proline=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[38],"Organ"))
  p_ME <- compare_means(Malate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME38 <- data.frame(Organ=p_ME$Organ,Malate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[39],"Organ"))
  p_ME <- compare_means(Citrate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME39 <- data.frame(Organ=p_ME$Organ,Citrate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[40],"Organ"))
  p_ME <- compare_means(UDP.glucuronic.acid~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME40 <- data.frame(Organ=p_ME$Organ,UDP.glucuronic.acid=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Group",predME[41],"Organ"))
  p_ME <- compare_means(Valine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME41 <- data.frame(Organ=p_ME$Organ,Valine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[42],"Organ"))
  p_ME <- compare_means(Cholesterol~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME42 <- data.frame(Organ=p_ME$Organ,Cholesterol=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[43],"Organ"))
  p_ME <- compare_means(Threonine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME43 <- data.frame(Organ=p_ME$Organ,Threonine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[44],"Organ"))
  p_ME <- compare_means(X3PD~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME44 <- data.frame(Organ=p_ME$Organ,X3PD=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[45],"Organ"))
  p_ME <- compare_means(dCMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME45 <- data.frame(Organ=p_ME$Organ,dCMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[46],"Organ"))
  p_ME <- compare_means(Lactate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME46 <- data.frame(Organ=p_ME$Organ,Lactate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[47],"Organ"))
  p_ME <- compare_means(Hypoxanthine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME47 <- data.frame(Organ=p_ME$Organ,Hypoxanthine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[48],"Organ"))
  p_ME <- compare_means(Glucose~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME48 <- data.frame(Organ=p_ME$Organ,Glucose=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[49],"Organ"))
  p_ME <- compare_means(Citrulline~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME49 <- data.frame(Organ=p_ME$Organ,Citrulline=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[50],"Organ"))
  p_ME <- compare_means(GABA~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME50 <- data.frame(Organ=p_ME$Organ,GABA=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Group",predME[51],"Organ"))
  p_ME <- compare_means(dTMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME51 <- data.frame(Organ=p_ME$Organ,dTMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[52],"Organ"))
  p_ME <- compare_means(dUMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME52 <- data.frame(Organ=p_ME$Organ,dUMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[53],"Organ"))
  p_ME <- compare_means(Xanthine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME53 <- data.frame(Organ=p_ME$Organ,Xanthine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[54],"Organ"))
  p_ME <- compare_means(Chondroitin~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME54 <- data.frame(Organ=p_ME$Organ,Chondroitin=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[55],"Organ"))
  p_ME <- compare_means(Isoleucine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME55 <- data.frame(Organ=p_ME$Organ,Isoleucine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[56],"Organ"))
  p_ME <- compare_means(X.E.E..Farnesyl.PP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME56 <- data.frame(Organ=p_ME$Organ,X.E.E..Farnesyl.PP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[57],"Organ"))
  p_ME <- compare_means(Deoxyadenosine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME57 <- data.frame(Organ=p_ME$Organ,Deoxyadenosine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[58],"Organ"))
  p_ME <- compare_means(XMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME58 <- data.frame(Organ=p_ME$Organ,XMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[59],"Organ"))
  p_ME <- compare_means(G6P~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME59 <- data.frame(Organ=p_ME$Organ,G6P=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[60],"Organ"))
  p_ME <- compare_means(dCDP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME60 <- data.frame(Organ=p_ME$Organ,dCDP=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Group",predME[61],"Organ"))
  p_ME <- compare_means(Spermine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME61 <- data.frame(Organ=p_ME$Organ,Spermine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[62],"Organ"))
  p_ME <- compare_means(Argininosuccinate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME62 <- data.frame(Organ=p_ME$Organ,Argininosuccinate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[63],"Organ"))
  p_ME <- compare_means(Dolichyl.phosphate.D.mannose~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME63 <- data.frame(Organ=p_ME$Organ,Dolichyl.phosphate.D.mannose=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[64],"Organ"))
  p_ME <- compare_means(AICAR~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME64 <- data.frame(Organ=p_ME$Organ,AICAR=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[65],"Organ"))
  p_ME <- compare_means(Fatty.Acid~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME65 <- data.frame(Organ=p_ME$Organ,Fatty.Acid=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[66],"Organ"))
  p_ME <- compare_means(Pyrimidine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME66 <- data.frame(Organ=p_ME$Organ,Pyrimidine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[67],"Organ"))
  p_ME <- compare_means(X.Glc.3..GlcNAc.2..Man.9..Asn.1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME67 <- data.frame(Organ=p_ME$Organ,X.Glc.3..GlcNAc.2..Man.9..Asn.1=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[68],"Organ"))
  p_ME <- compare_means(X.GlcNAc.4..Man.3..Asn.1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME68 <- data.frame(Organ=p_ME$Organ,X.GlcNAc.4..Man.3..Asn.1=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[69],"Organ"))
  p_ME <- compare_means(X.Gal.2..GlcA.1..Xyl.1..Ser.1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  p_ME69 <- data.frame(Organ=p_ME$Organ,X.Gal.2..GlcA.1..Xyl.1..Ser.1=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[70],"Organ"))
  p_ME <- compare_means(X.Gal.1..GlcNAc.1..Man.1..Ser.Thr.1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
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
    data <- FetchData(scRNA, vars = c("Group",predME[1],"Organ"))
    p_ME <- compare_means(AMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME1 <- data.frame(Organ=p_ME$Organ,AMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[2],"Organ"))
    p_ME <- compare_means(Pyruvate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME2 <- data.frame(Organ=p_ME$Organ,Pyruvate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[3],"Organ"))
    p_ME <- compare_means(Acetyl.CoA~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME3 <- data.frame(Organ=p_ME$Organ,Acetyl.CoA=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[4],"Organ"))
    p_ME <- compare_means(Glutamate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME4 <- data.frame(Organ=p_ME$Organ,Glutamate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[5],"Organ"))
    p_ME <- compare_means(X2OG~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME5 <- data.frame(Organ=p_ME$Organ,X2OG=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[6],"Organ"))
    p_ME <- compare_means(Oxaloacetate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME6 <- data.frame(Organ=p_ME$Organ,Oxaloacetate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[7],"Organ"))
    p_ME <- compare_means(Glycine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME7 <- data.frame(Organ=p_ME$Organ,Glycine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[8],"Organ"))
    p_ME <- compare_means(Succinate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME8 <- data.frame(Organ=p_ME$Organ,Succinate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[9],"Organ"))
    p_ME <- compare_means(UDP.N.acetylglucosamine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME9 <- data.frame(Organ=p_ME$Organ,UDP.N.acetylglucosamine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[10],"Organ"))
    p_ME <- compare_means(lysine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME10 <- data.frame(Organ=p_ME$Organ,lysine=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[11],"Organ"))
    p_ME <- compare_means(Aspartate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME11 <- data.frame(Organ=p_ME$Organ,Aspartate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[12],"Organ"))
    p_ME <- compare_means(Glutathione~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME12 <- data.frame(Organ=p_ME$Organ,Glutathione=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[13],"Organ"))
    p_ME <- compare_means(Arginine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME13 <- data.frame(Organ=p_ME$Organ,Arginine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[14],"Organ"))
    p_ME <- compare_means(Glutamine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME14 <- data.frame(Organ=p_ME$Organ,Glutamine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[15],"Organ"))
    p_ME <- compare_means(Serine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME15 <- data.frame(Organ=p_ME$Organ,Serine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[16],"Organ"))
    p_ME <- compare_means(Methionine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME16 <- data.frame(Organ=p_ME$Organ,Methionine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[17],"Organ"))
    p_ME <- compare_means(Ornithine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME17 <- data.frame(Organ=p_ME$Organ,Ornithine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[18],"Organ"))
    p_ME <- compare_means(Phenylalanine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME18 <- data.frame(Organ=p_ME$Organ,Phenylalanine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[19],"Organ"))
    p_ME <- compare_means(Tyrosine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME19 <- data.frame(Organ=p_ME$Organ,Tyrosine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[20],"Organ"))
    p_ME <- compare_means(Succinyl.CoA~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME20 <- data.frame(Organ=p_ME$Organ,Succinyl.CoA=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[21],"Organ"))
    p_ME <- compare_means(Cysteine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME21 <- data.frame(Organ=p_ME$Organ,Cysteine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[22],"Organ"))
    p_ME <- compare_means(B.Alanine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME22 <- data.frame(Organ=p_ME$Organ,B.Alanine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[23],"Organ"))
    p_ME <- compare_means(Propanoyl.CoA~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME23 <- data.frame(Organ=p_ME$Organ,Propanoyl.CoA=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[24],"Organ"))
    p_ME <- compare_means(Glucose.1.phosphate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME24 <- data.frame(Organ=p_ME$Organ,Glucose.1.phosphate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[25],"Organ"))
    p_ME <- compare_means(UMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME25 <- data.frame(Organ=p_ME$Organ,UMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[26],"Organ"))
    p_ME <- compare_means(Dolichyl.phosphate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME26 <- data.frame(Organ=p_ME$Organ,Dolichyl.phosphate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[27],"Organ"))
    p_ME <- compare_means(CDP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME27 <- data.frame(Organ=p_ME$Organ,CDP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[28],"Organ"))
    p_ME <- compare_means(Choline~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME28 <- data.frame(Organ=p_ME$Organ,Choline=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[29],"Organ"))
    p_ME <- compare_means(G3P~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME29 <- data.frame(Organ=p_ME$Organ,G3P=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[30],"Organ"))
    p_ME <- compare_means(PRPP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME30 <- data.frame(Organ=p_ME$Organ,PRPP=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[31],"Organ"))
    p_ME <- compare_means(Fumarate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME31 <- data.frame(Organ=p_ME$Organ,Fumarate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[32],"Organ"))
    p_ME <- compare_means(Leucine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME32 <- data.frame(Organ=p_ME$Organ,Leucine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[33],"Organ"))
    p_ME <- compare_means(IMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME33 <- data.frame(Organ=p_ME$Organ,IMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[34],"Organ"))
    p_ME <- compare_means(Putrescine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME34 <- data.frame(Organ=p_ME$Organ,Putrescine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[35],"Organ"))
    p_ME <- compare_means(Histidine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME35 <- data.frame(Organ=p_ME$Organ,Histidine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[36],"Organ"))
    p_ME <- compare_means(GMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME36 <- data.frame(Organ=p_ME$Organ,GMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[37],"Organ"))
    p_ME <- compare_means(Proline~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME37 <- data.frame(Organ=p_ME$Organ,Proline=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[38],"Organ"))
    p_ME <- compare_means(Malate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME38 <- data.frame(Organ=p_ME$Organ,Malate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[39],"Organ"))
    p_ME <- compare_means(Citrate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME39 <- data.frame(Organ=p_ME$Organ,Citrate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[40],"Organ"))
    p_ME <- compare_means(UDP.glucuronic.acid~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME40 <- data.frame(Organ=p_ME$Organ,UDP.glucuronic.acid=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[41],"Organ"))
    p_ME <- compare_means(Valine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME41 <- data.frame(Organ=p_ME$Organ,Valine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[42],"Organ"))
    p_ME <- compare_means(Cholesterol~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME42 <- data.frame(Organ=p_ME$Organ,Cholesterol=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[43],"Organ"))
    p_ME <- compare_means(Threonine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME43 <- data.frame(Organ=p_ME$Organ,Threonine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[44],"Organ"))
    p_ME <- compare_means(X3PD~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME44 <- data.frame(Organ=p_ME$Organ,X3PD=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[45],"Organ"))
    p_ME <- compare_means(dCMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME45 <- data.frame(Organ=p_ME$Organ,dCMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[46],"Organ"))
    p_ME <- compare_means(Lactate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME46 <- data.frame(Organ=p_ME$Organ,Lactate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[47],"Organ"))
    p_ME <- compare_means(Hypoxanthine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME47 <- data.frame(Organ=p_ME$Organ,Hypoxanthine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[48],"Organ"))
    p_ME <- compare_means(Glucose~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME48 <- data.frame(Organ=p_ME$Organ,Glucose=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[49],"Organ"))
    p_ME <- compare_means(Citrulline~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME49 <- data.frame(Organ=p_ME$Organ,Citrulline=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[50],"Organ"))
    p_ME <- compare_means(GABA~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME50 <- data.frame(Organ=p_ME$Organ,GABA=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[51],"Organ"))
    p_ME <- compare_means(dTMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME51 <- data.frame(Organ=p_ME$Organ,dTMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[52],"Organ"))
    p_ME <- compare_means(dUMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME52 <- data.frame(Organ=p_ME$Organ,dUMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[53],"Organ"))
    p_ME <- compare_means(Xanthine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME53 <- data.frame(Organ=p_ME$Organ,Xanthine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[54],"Organ"))
    p_ME <- compare_means(Chondroitin~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME54 <- data.frame(Organ=p_ME$Organ,Chondroitin=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[55],"Organ"))
    p_ME <- compare_means(Isoleucine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME55 <- data.frame(Organ=p_ME$Organ,Isoleucine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[56],"Organ"))
    p_ME <- compare_means(X.E.E..Farnesyl.PP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME56 <- data.frame(Organ=p_ME$Organ,X.E.E..Farnesyl.PP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[57],"Organ"))
    p_ME <- compare_means(Deoxyadenosine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME57 <- data.frame(Organ=p_ME$Organ,Deoxyadenosine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[58],"Organ"))
    p_ME <- compare_means(XMP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME58 <- data.frame(Organ=p_ME$Organ,XMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[59],"Organ"))
    p_ME <- compare_means(G6P~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME59 <- data.frame(Organ=p_ME$Organ,G6P=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[60],"Organ"))
    p_ME <- compare_means(dCDP~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME60 <- data.frame(Organ=p_ME$Organ,dCDP=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[61],"Organ"))
    p_ME <- compare_means(Spermine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME61 <- data.frame(Organ=p_ME$Organ,Spermine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[62],"Organ"))
    p_ME <- compare_means(Argininosuccinate~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME62 <- data.frame(Organ=p_ME$Organ,Argininosuccinate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[63],"Organ"))
    p_ME <- compare_means(Dolichyl.phosphate.D.mannose~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME63 <- data.frame(Organ=p_ME$Organ,Dolichyl.phosphate.D.mannose=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[64],"Organ"))
    p_ME <- compare_means(AICAR~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME64 <- data.frame(Organ=p_ME$Organ,AICAR=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[65],"Organ"))
    p_ME <- compare_means(Fatty.Acid~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME65 <- data.frame(Organ=p_ME$Organ,Fatty.Acid=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[66],"Organ"))
    p_ME <- compare_means(Pyrimidine~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME66 <- data.frame(Organ=p_ME$Organ,Pyrimidine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[67],"Organ"))
    p_ME <- compare_means(X.Glc.3..GlcNAc.2..Man.9..Asn.1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME67 <- data.frame(Organ=p_ME$Organ,X.Glc.3..GlcNAc.2..Man.9..Asn.1=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[68],"Organ"))
    p_ME <- compare_means(X.GlcNAc.4..Man.3..Asn.1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME68 <- data.frame(Organ=p_ME$Organ,X.GlcNAc.4..Man.3..Asn.1=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[69],"Organ"))
    p_ME <- compare_means(X.Gal.2..GlcA.1..Xyl.1..Ser.1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    p_ME69 <- data.frame(Organ=p_ME$Organ,X.Gal.2..GlcA.1..Xyl.1..Ser.1=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[70],"Organ"))
    p_ME <- compare_means(X.Gal.1..GlcNAc.1..Man.1..Ser.Thr.1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
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

r2 <-as.matrix(read.csv('logp_score_all.csv', header = T, row.names = 1))
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








