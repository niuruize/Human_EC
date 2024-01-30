################################################################################
###---Human_Chi_mac_EC_integration
################################################################################
library(Seurat);library(tidyverse);library(patchwork);library(ggplot2);library(DoubletFinder)
counts <- Matrix::readMM('nuclei_count_matrix.mtx.gz')
counts <- counts[,-1]
meta <-read.csv('metadata_nuclei.tsv', header = TRUE,sep = "\t",row.names = 1)
gene <-read.csv('genes.tsv', header = F,sep = "\t")
colnames(counts) <- rownames(meta)
rownames(counts) <- gene$V1
E_MTAB_8230 <- CreateSeuratObject(counts = counts, meta.data = meta,min.cells = 3, min.features = 200) 
E_MTAB_8230 <- AddMetaData(E_MTAB_8230, metadata = meta$Sample_ID,col.name = "SampleID")
E_MTAB_8230 <- AddMetaData(E_MTAB_8230, metadata = meta$Cell_class,col.name = "celltype")
E_MTAB_8230$Age <- plyr::mapvalues(x = as.character(E_MTAB_8230$SampleID), from = current.cluster.ids, to = new.cluster.ids)
E_MTAB_8230[["percent.mt"]] <- E_MTAB_8230$Percentage_mito
summary(E_MTAB_8230[[]]$percent.mt)
E_MTAB_8230[["percent.rb"]] <- PercentageFeatureSet(E_MTAB_8230,pattern = "^RP[SL]")
summary(E_MTAB_8230[[]]$percent.rb)
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(E_MTAB_8230))
E_MTAB_8230[["percent.HB"]] <- PercentageFeatureSet(E_MTAB_8230,features = HB.genes)
summary(E_MTAB_8230[[]]$percent.HB)
summary(E_MTAB_8230[[]]$nFeature_RNA)
summary(E_MTAB_8230[[]]$nCount_RNA)
VlnPlot(E_MTAB_8230, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
#E_MTAB_8230 <- subset(E_MTAB_8230, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 &  percent.mt < 10 & percent.HB < 1)
#VlnPlot(E_MTAB_8230, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
#plotE_MTAB_8230_1 <- FeatureScatter(E_MTAB_8230, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "SampleID")
#plotE_MTAB_8230_2 <- FeatureScatter(E_MTAB_8230, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "SampleID")
#plotE_MTAB_8230_1 + plotE_MTAB_8230_2
E_MTAB_8230[["datasets"]] <- "E_MTAB_8230_Human_Chimp_Macaque "
save(E_MTAB_8230 ,file="E_MTAB_8230.RData")
###---integration
load("~/aging/EC/E_MTAB_8230.RData")
table(E_MTAB_8230$Sample_ID)
scRNA=E_MTAB_8230
scRNA <- SCTransform(scRNA)
###--PCA
scRNA <- RunPCA(scRNA, npcs=30, verbose=FALSE)
ElbowPlot(scRNA, ndims = 30)
###
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=20)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:20) %>% RunUMAP(reduction="harmony", dims = 1:20)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:20) %>% FindClusters(dims = 1:20, resolution = 0.8)
save(scRNA,file="Hu_Chi_Mac.RData")







