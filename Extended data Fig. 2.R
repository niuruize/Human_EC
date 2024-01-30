################################################################################
# Extended data Fig. 2A
################################################################################
load('Human_EC.RData')
table(scRNA$Organ)  
av<-AverageExpression(scRNA,group.by = "Region", assays = "SCT")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
write.csv(cor(av[cg,],method = "spearman"),"cor_Region.csv") #保存结果
pdf("cor_Region.pdf", width = 8,height = 7.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()

################################################################################
# Extended data Fig. 2B
################################################################################
library(dendextend);library(circlize);library(Seurat)
load('Human_EC.RData')
##
cluster.averages <- AverageExpression(scRNA,group.by = "Region")
B.exp <- cluster.averages[["SCT"]]
hc <- as.dendrogram(hclust(dist(t(B.exp)), "ave"))
hc <- hc %>% color_branches(k = 20) %>% color_labels(k = 38)
##
pdf("cluster_region.pdf", width = 10,height = 3)
plot(hc)
dev.off()
##
pdf("circ_cluster_region.pdf", width = 8,height = 8)
circlize_dendrogram(hc,labels_track_height = NA,dend_track_height = 0.3)
dev.off()

################################################################################
# Extended data Fig. 2C
################################################################################
library(ggplot2);library(scRNAtoolVis)
markers <- read_xlsx("annogene.xlsx")
p1 = jjDotPlot(scRNA, gene = annoGene$gene,id = 'Organ',xtree = T,ytree = T,rescale = T,rescale.min = 0,rescale.max = 1)
ggsave(filename = "EC_marker_heatmap.pdf", plot = p1, device = 'pdf', width = 21, height = 21, units = 'cm')

################################################################################
# Extended data Fig. 2D
################################################################################
# GSEA
library(org.Hs.eg.db);library(clusterProfiler);library(enrichplot);library(ggplot2);library(enrichplot)
markers = read.table("markers.txt",header = T) #读入txt
markers = subset(markers, markers$p_val_adj < 0.05)
celltype = c("Aorta","Bladder","Brain","Decidua","Esophagus","Heart","Intestines","Kidney","Liver","Lung","MG","Muscle","Ovary","Pancreas","PIT","Placenta","Retina","SAT","SC","Skin","Testis","Thyroid","Trachea","Uterus","VAT")
GSEA_organ <- data.frame(GO = numeric(0))
for (i in 1:length(celltype)) {
  df <- subset(markers, markers$cluster==celltype[i])
  df_id <- bitr(df$SYMBOL, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
  df_all <- merge(df,df_id,by="SYMBOL",all=F)
  #GSEA
  df_all_sort <- df_all[order(df_all$avg_log2FC, decreasing = T),]
  gene_fc = df_all_sort$avg_log2FC
  names(gene_fc) <- df_all_sort$ENTREZID
  #GO富集
  GO <- gseGO(gene_fc, ont = "BP",OrgDb = org.Hs.eg.db,keyType = "ENTREZID",pvalueCutoff = 0.05,pAdjustMethod = "BH")
  #数据导出
  sortGO <- GO[order(GO$NES, decreasing = T),]
  sortGO <- data.frame(sortGO$Description, sortGO$NES)
  names(sortGO) <- c("GO", celltype[i])
  GSEA_organ <- merge(GSEA_organ, sortGO,by="GO",all=TRUE)
}
write.csv(GSEA_organ,"GSEA_organ.csv", row.names = F)







