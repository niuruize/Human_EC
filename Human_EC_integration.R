################################################################################
###---Human_EC_integration
################################################################################
###
###--------------------------------------harmony integration---------------------------------------####
#library packages
library(Seurat);library(harmony);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

#load data
load("~/Downloads/scRNA/EC/adipose/Adipose_Endo1.RData")
load("~/Downloads/scRNA/EC/aorta/Aorta_Endo1.RData")
load("~/Downloads/scRNA/EC/Bladder/Bladder_Endo1.RData")
Brain_Endo$Region <- as.factor(as.character(Brain_Endo$Region))
load("~/Downloads/scRNA/EC/brain/Brain_Endo1.RData")
load("~/Downloads/scRNA/EC/Esophagus/Esophagus_Endo1.RData")
load("~/Downloads/scRNA/EC/gland/MG_Endo1.RData")
load("~/Downloads/scRNA/EC/heart/Heart_Endo1.RData")
load("~/Downloads/scRNA/EC/Intestinal/GI_Endo1.RData")
load("~/Downloads/scRNA/EC/kidney/Kidney_Endo1.RData")
load("~/Downloads/scRNA/EC/liver/Liver_Endo1.RData")
load("~/Downloads/scRNA/EC/lung/Lung_Endo1.RData")
load("~/Downloads/scRNA/EC/maternal_fetal/decidua/Decidua_Endo1.RData")
load("~/Downloads/scRNA/EC/maternal_fetal/placenta/Placenta_Endo1.RData")
load("~/Downloads/scRNA/EC/muscle/Muscle_Endo1.RData")
load("~/Downloads/scRNA/EC/ovary/Ovary_Endo1.RData")
load("~/Downloads/scRNA/EC/pancreas/Pancreas_Endo1.RData")
load("~/Downloads/scRNA/EC/PIT/PIT_Endo1.RData")
load("~/Downloads/scRNA/EC/SC/SC_Endo1.RData")
load("~/Downloads/scRNA/EC/skin/Skin_Endo1.RData")
load("~/Downloads/scRNA/EC/testis/Testis_Endo1.RData")
load("~/Downloads/scRNA/EC/thyroid/Thyroid_Endo1.RData")
load("~/Downloads/scRNA/EC/Trachea/Trachea_Endo1.RData")
load("~/Downloads/scRNA/EC/uterus/Uterus_Endo1.RData")
load("~/Downloads/scRNA/EC/retina/retina_Endo1.RData")

#merge data
scRNA = merge(Adipose_Endo, y=c(Aorta_Endo,Bladder_Endo,Brain_Endo,Decidua_Endo,Esophagus_Endo,GI_Endo,Heart_Endo,Kidney_Endo,Liver_Endo,Lung_Endo,
                                MG_Endo,Muscle_Endo,Ovary_Endo,Pancreas_Endo,PIT_Endo,Placenta_Endo,Retina_Endo,SC_Endo,Skin_Endo,Testis_Endo,Thyroid_Endo,
                                Trachea_Endo,Uterus_Endo))
table(scRNA$Organ)
table(scRNA$Region)
table(scRNA$Organ,scRNA$Region)
save(scRNA,file="Endo_merge.RData")

###
scRNA <- SCTransform(scRNA)
###--PCA
scRNA <- RunPCA(scRNA, npcs=30, verbose=FALSE)
ElbowPlot(scRNA, ndims = 30)
##
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=20)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:20) %>% RunUMAP(reduction="harmony", dims = 1:20)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:20) %>% FindClusters(dims = 1:20, resolution = 0.4)
save(scRNA,file="Human_EC.RData")

DimPlot(scRNA,group.by = "SampleID",reduction = "umap")
DimPlot(scRNA,group.by = "seurat_clusters",reduction = "umap",label = TRUE)
DimPlot(scRNA,group.by = "SCT_snn_res.1.2",reduction = "umap",label = TRUE)
DimPlot(scRNA,group.by = "SampleID", reduction = "umap", split.by = "Region")

table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["Sex"]])
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["celltype"]])
table(scRNA@meta.data[["Age"]],scRNA@meta.data[["seurat_clusters"]])
table(scRNA@meta.data[["seurat_clusters"]],scRNA@meta.data[["Age"]])

#
DefaultAssay(scRNA) <- "SCT"
Idents(scRNA) = "Organ"
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = FALSE, test.use = "wilcox")  ##耗时久
write.table(markers,file="markers.txt",quote=F,sep="\t",row.names=F,col.names=T)
markers_df = markers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_log2FC)
write.table(markers_df,file="markers500.txt",quote=F,sep="\t",row.names=F,col.names=T)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10,file="markers10.txt",quote=F,sep="\t",row.names=F,col.names=T)
Heatmap <- DoHeatmap(scRNA,features = top10$gene) 
ggsave(filename = "Heatmap10.pdf", plot = Heatmap, device = 'pdf', width = 60, height = 60, units = 'cm')
top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
subobj <- subset(scRNA, downsample = 500)
options(repr.plot.width = 13, repr.plot.height=13)
Heatmap <- DoHeatmap(subobj,features = top5$gene,disp.min=-2.5, disp.max=2.5,group.colors = c("#00BFC4","#AB82FF","#00CD00","#C77CFF"))+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))  
ggsave(filename = "Heatmap5_1.pdf", plot = Heatmap, device = 'pdf', width = 60, height = 60, units = 'cm')
Heatmap <- DoHeatmap(scRNA,features = top5$gene,disp.min=-2.5, disp.max=2.5, group.colors = c("#00BFC4","#AB82FF","#00CD00","#C77CFF"))+
  scale_fill_gradientn(colors = c("white","grey","firebrick3")) 
ggsave(filename = "Heatmap5_2.pdf", plot = Heatmap, device = 'pdf', width = 60, height = 60, units = 'cm')

# correlation analysis
table(scRNA$celltype)  
av<-AverageExpression(scRNA,group.by = "celltype", assays = "RNA")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_celltype.pdf", width = 4,height = 3.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_celltype.csv") #保存结果
#
table(scRNA$celltype1)  
av<-AverageExpression(scRNA,group.by = "celltype1", assays = "RNA")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_celltype1.pdf", width = 7,height = 6.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_celltype1.csv") #保存结果
#
table(scRNA$Organ)  
av<-AverageExpression(scRNA,group.by = "Organ", assays = "SCT")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
write.csv(cor(av[cg,],method = "spearman"),"cor_Organ.csv") #保存结果
pdf("cor_Organ.pdf", width = 8,height = 7.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()
#
table(scRNA$Organ)  
av<-AverageExpression(scRNA,group.by = "Region", assays = "SCT")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
write.csv(cor(av[cg,],method = "spearman"),"cor_Region.csv") #保存结果
pdf("cor_Region.pdf", width = 8,height = 7.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()

##
DefaultAssay(scRNA)="RNA"
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
save(scRNA,file="Human_EC.RData")

######------大型数据保存之前最好清空scale.data和data
### Before saving
#scRNA_clean <- Seurat_object ## 10W cells, 15G
scRNA_clean <- scRNA
scRNA_clean@assays$RNA@data <- as.matrix(0)
scRNA_clean@assays$RNA@scale.data <- as.matrix(0)
save(scRNA_clean, file = "Seurat.RData") ### ~0.2-0.3G

###loading
load("Seurat.RData")
DefaultAssay(scRNA_clean) <- "RNA"
scRNA_clean@assays$RNA@data <- scRNA_clean@assays$RNA@counts
scRNA_clean <- NormalizeData(scRNA_clean)





