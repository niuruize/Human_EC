################################################################################
###---Mouse_EC_integration
################################################################################
###
library(Seurat);library(tidyverse);library(patchwork);library(ggplot2);library(DoubletFinder)
#if (!requireNamespace("BiocManager", quietly = TRUE)) 
#  install.packages("BiocManager") 
#BiocManager::install("Nebulosa") 
counts <- read.table("raw_count_matrix_brain.txt", header = TRUE, sep = "\t",row.names = 1)
brain <- CreateSeuratObject(counts = counts, project = "brain",min.cells = 3, min.features = 200) 
brain$Organ <- "brain"
counts <- read.table("raw_count_matrix_colon.txt", header = TRUE, sep = "\t",row.names = 1)
colon <- CreateSeuratObject(counts = counts, project = "colon",min.cells = 3, min.features = 200) 
colon$Organ <- "Colon"
counts <- read.table("raw_count_matrix_heart.txt", header = TRUE, sep = "\t",row.names = 1)
heart <- CreateSeuratObject(counts = counts, project = "heart",min.cells = 3, min.features = 200) 
heart$Organ <- "heart"
counts <- read.table("raw_count_matrix_kidney.txt", header = TRUE, sep = "\t",row.names = 1)
kidney <- CreateSeuratObject(counts = counts, project = "kidney",min.cells = 3, min.features = 200) 
kidney$Organ <- "kidney"
counts <- read.table("raw_count_matrix_liver.txt", header = TRUE, sep = "\t",row.names = 1)
liver <- CreateSeuratObject(counts = counts, project = "Liver",min.cells = 3, min.features = 200) 
liver$Organ <- "Liver"
counts <- read.table("raw_count_matrix_lung.txt", header = TRUE, sep = "\t",row.names = 1)
lung <- CreateSeuratObject(counts = counts, project = "lung",min.cells = 3, min.features = 200) 
lung$Organ <- "lung"
counts <- read.table("raw_count_matrix_muscle_EDL.txt", header = TRUE, sep = "\t",row.names = 1)
EDL <- CreateSeuratObject(counts = counts, project = "EDL",min.cells = 3, min.features = 200) 
EDL$Organ <- "EDL"
counts <- read.table("raw_count_matrix_muscle_soleus.txt", header = TRUE, sep = "\t",row.names = 1)
soleus <- CreateSeuratObject(counts = counts, project = "soleus",min.cells = 3, min.features = 200) 
soleus$Organ <- "soleus"
counts <- read.table("raw_count_matrix_small_intestine.txt", header = TRUE, sep = "\t",row.names = 1)
small_intestine <- CreateSeuratObject(counts = counts, project = "small_intestine",min.cells = 3, min.features = 200) 
small_intestine$Organ <- "small_intestine"
counts <- read.table("raw_count_matrix_spleen.txt", header = TRUE, sep = "\t",row.names = 1)
spleen <- CreateSeuratObject(counts = counts, project = "spleen",min.cells = 3, min.features = 200) 
spleen$Organ <- "spleen"
counts <- read.table("raw_count_matrix_testis.txt", header = TRUE, sep = "\t",row.names = 1)
testis <- CreateSeuratObject(counts = counts, project = "testis",min.cells = 3, min.features = 200) 
testis$Organ <- "testis"
# merge
E_MTAB_8077 = merge(brain,y=c(colon,EDL,heart,kidney,liver,lung,small_intestine,soleus,spleen,testis))
E_MTAB_8077[["percent.mt"]] <- PercentageFeatureSet(E_MTAB_8077,pattern = "^mt-")
summary(E_MTAB_8077[[]]$percent.mt)
summary(E_MTAB_8077[[]]$nFeature_RNA)
summary(E_MTAB_8077[[]]$nCount_RNA)
VlnPlot(E_MTAB_8077, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E_MTAB_8077 <- subset(E_MTAB_8077, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 &  percent.mt < 5)
E_MTAB_8077[["datasets"]] <- "E_MTAB_8077_mouse"
E_MTAB_8077 <- SCTransform(E_MTAB_8077, vars.to.regress = "percent.mt", verbose = FALSE)
E_MTAB_8077 <- RunPCA(E_MTAB_8077)
E_MTAB_8077 <- FindNeighbors(E_MTAB_8077, dims = 1:30)
E_MTAB_8077 <- FindClusters(E_MTAB_8077, resolution = 0.8)
E_MTAB_8077 <- RunUMAP(E_MTAB_8077, dims = 1:30)
sweep.res.list <- paramSweep_v3(E_MTAB_8077, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(E_MTAB_8077)*8*1e-6
homotypic.prop <- modelHomotypic(E_MTAB_8077$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(E_MTAB_8077))) ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
E_MTAB_8077 <- doubletFinder_v3(E_MTAB_8077, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(E_MTAB_8077@meta.data)[ncol(E_MTAB_8077@meta.data)]="DoubletFinder"
DimPlot(E_MTAB_8077,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
rm('plotE_MTAB_8077_1','plotE_MTAB_8077_2','HB.genes','counts','sweep.res.list','sweep.stats','DoubletRate','nExp_poi','nExp_poi.agj','pk_good','pk_v','homotypic.prop','bcmvn')
cellinfo <- subset(E_MTAB_8077@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Organ","DoubletFinder","datasets"))
E_MTAB_8077 <- CreateSeuratObject(E_MTAB_8077@assays$RNA@counts, meta.data = cellinfo)
save(E_MTAB_8077 ,file="E_MTAB_8077.RData")
# integration
scRNA = subset(E_MTAB_8077, DoubletFinder == "Singlet")
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=30, verbose=FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- RunHarmony(scRNA, group.by.vars="Organ",assay.use="SCT",max.iter.harmony=20)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:20) %>% RunUMAP(reduction="harmony", dims = 1:20)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:20) %>% FindClusters(dims = 1:20, resolution = 0.2)
save(scRNA,file="Mouse_EC_seurat.RData")
# QC
scRNA = scRNA[,scRNA$seurat_clusters %in% c("0","1","2","3","4","5","6","7","8","9","10")]
cellinfo <- subset(scRNA@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Organ","DoubletFinder","datasets"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=30, verbose=FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- RunHarmony(scRNA, group.by.vars="Organ",assay.use="SCT",max.iter.harmony=20)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:20) %>% RunUMAP(reduction="harmony", dims = 1:20)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:20) %>% FindClusters(dims = 1:20, resolution = 0.2)
save(scRNA,file="Mouse_EC.RData")








