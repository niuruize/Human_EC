################################################################################
# Figure 6 scDRS
################################################################################
# extract expression counts matrix
library(Seurat);library(tidyverse);library(cowplot);library(patchwork);library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db);library(MySeuratWrappers);library(scRNAtoolVis);library(readxl);library(harmony)
load("Human_EC.RData")
scRNA$celltype <- as.character(scRNA$celltype)
scRNA$celltype_label <- as.character(scRNA$celltype_label)
scRNA$celltype_label2 <- as.character(scRNA$celltype_label2)
scRNA$Organ <- as.character(scRNA$Organ)
scRNA$Region <- as.character(scRNA$Region)
count <- scRNA@assays$RNA@counts
#scRNA$Group_celltype <- paste(scRNA$Group, scRNA$celltype, sep = "_")
scRNA2<-CreateSeuratObject(counts=count, meta.data = scRNA@meta.data)
SaveH5Seurat(scRNA2,filename="EC.h5seurat",overwrite = T)
Convert("EC.h5seurat",dest="h5ad",assay="RNA",overwrite = T)

EC_cov <- subset(scRNA@meta.data, select= c("Sex","Age","percent.mt","percent.rb","nCount_RNA"))
EC_cov$Sex = factor(EC_cov$Sex)
EC_cov$Age = as.numeric(EC_cov$Age)
EC_cov$percent.mt = as.numeric(EC_cov$percent.mt)
EC_cov$percent.rb = as.numeric(EC_cov$percent.rb)
EC_cov$nCount_RNA = as.numeric(EC_cov$nCount_RNA)
write.table(EC_cov, file='EC_cov.tsv',sep = '\t')

# save metadata table:
scRNA$barcode <- colnames(scRNA)
scRNA$UMAP_1 <- scRNA@reductions$umap@cell.embeddings[,1]
scRNA$UMAP_2 <- scRNA@reductions$umap@cell.embeddings[,2]
write.csv(scRNA@meta.data, file='metadata.csv', quote=F, row.names=F)

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(scRNA@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(scRNA, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.mtx')

# write gene names
write.table(data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',quote=F,row.names=F,col.names=F)

# scDRS analysis used python
nohup scdrs compute-score \
    --h5ad-file EC.h5ad \
    --h5ad-species human \
    --gs-file Lung/processed_geneset.gs \
    --cov-file EC_cov.tsv \
    --gs-species human \
    --flag-filter-data True \
    --flag-raw-count True \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True \
    --out-folder EC/ &

#=====================================================================
#                       Figure 6 scDRS visualization
#=====================================================================
# Figure 6A-F
# The code for the visualization is shown in scDRS.ipynb





