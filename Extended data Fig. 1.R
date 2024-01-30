#loading R packages
library(Seurat);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

##########################################################################################################################
# Adipose
##########################################################################################################################
meta <-read.csv('GSE176171_cell_metadata.tsv.gz', header = TRUE,sep = "\t")
meta <- meta[1:137684,]
GSE176171 <- Read10X(data.dir = "/Adipose/PMID35296864/GSE176171")
GSE176171 <- CreateSeuratObject(counts = GSE176171, project = "GSE176171", min.cells = 3, min.features = 200) 
GSE176171 <- AddMetaData(GSE176171, metadata = meta$biosample_id,col.name = "SampleID")
GSE176171 <- AddMetaData(GSE176171, metadata = meta$sex__ontology_label,col.name = "Sex")
GSE176171 <- AddMetaData(GSE176171, metadata = meta$cell_type__custom,col.name = "celltype")
GSE176171 <- AddMetaData(GSE176171, metadata = meta$bmi,col.name = "BMI")
GSE176171 <- AddMetaData(GSE176171, metadata = meta$bmi__group,col.name = "BMI_group")
GSE176171 <- AddMetaData(GSE176171, metadata = meta$fat__type,col.name = "Region")
GSE176171 <- AddMetaData(GSE176171, metadata = meta$organism_age,col.name = "Age")
GSE176171[["percent.mt"]] <- PercentageFeatureSet(GSE176171,pattern = "^MT-")
summary(GSE176171[[]]$percent.mt)
GSE176171[["percent.rb"]] <- PercentageFeatureSet(GSE176171,pattern = "^RP[SL]")
summary(GSE176171[[]]$percent.rb)
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(GSE176171))
GSE176171[["percent.HB"]] <- PercentageFeatureSet(GSE176171,features = HB.genes)
summary(GSE176171[[]]$percent.HB)
summary(GSE176171[[]]$nFeature_RNA)
summary(GSE176171[[]]$nCount_RNA)
VlnPlot(GSE176171, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
current.cluster.ids <- c("female","male")
new.cluster.ids <- c("F","M") 
GSE176171$Sex <- plyr::mapvalues(x = as.character(GSE176171$Sex), from = current.cluster.ids, to = new.cluster.ids)
current.cluster.ids <- c("24","29","35","36","41","51","53","57","58","68","73")
new.cluster.ids <- c("Adult","Adult","Adult","Adult","Adult","Adult","Adult","Adult","Adult","Aging","Aging") 
GSE176171$Group <- plyr::mapvalues(x = as.character(GSE176171$Age), from = current.cluster.ids, to = new.cluster.ids)
save(GSE176171 ,file="GSE176171.RData")

# loading data
load('GSE176171.RData')

# integrated data
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
save(scRNA,file="Adipose_seurat.RData")
#
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE)  ##耗时久
head(Idents(scRNA), 5)
save(scRNA,file="Adipose_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","DCN","COL1A2","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)
#
table(scRNA$seurat_clusters)
Adipose_Endo <- scRNA[,scRNA$seurat_clusters %in% c("9","11","12","17","26")]
Adipose_Endo <- CreateSeuratObject(Adipose_Endo@assays$RNA@counts, meta.data = Adipose_Endo@meta.data)
Adipose_Endo[["Organ"]] <- "Adipose"
save(Adipose_Endo,file="Adipose_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)

##########################################################################################################################
# Aorta
##########################################################################################################################
# loading data
load("~/Downloads/scRNA/EC/aorta/GSE155468/GSE155468_Aorta.RData")
load("~/Downloads/scRNA/EC/aorta/GSE166676/GSE166676_Aorta.RData")
load("~/Downloads/scRNA/EC/aorta/GSE207784/GSE207784_Aorta.RData")
GSE155468_Aorta = subset(GSE155468_Aorta, DoubletFinder == "Singlet")
GSE166676_Aorta = subset(GSE166676_Aorta, DoubletFinder == "Singlet")
GSE207784_Aorta = subset(GSE207784_Aorta, DoubletFinder == "Singlet")

# integrated data
scRNA = merge(GSE155468_Aorta, y=c(GSE166676_Aorta,GSE207784_Aorta))
rm(GSE155468_Aorta,GSE166676_Aorta,GSE207784_Aorta)
save(scRNA,file="Aorta_merge.RData")
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:30)
scRNA <- RunTSNE(scRNA, dims = 1:30,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="Aorta_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","DCN","COL1A2","COL1A1","ACTA2","PECAM1","VWF","CDH5","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)
table(scRNA$seurat_clusters)
Aorta_Endo <- scRNA[,scRNA$seurat_clusters %in% c("10","12","20")]
Aorta_Endo <- CreateSeuratObject(Aorta_Endo@assays$RNA@counts, meta.data = Aorta_Endo@meta.data)
Aorta_Endo[["Region"]] <- "Aorta"
Aorta_Endo[["Organ"]] <- "Aorta"
save(Aorta_Endo,file="Aorta_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Bladder
##########################################################################################################################
# loading data
load("~/Downloads/scRNA/EC/Bladder/GSE129845/GSE129845_Bladder.RData")
load("~/Downloads/scRNA/EC/Bladder/GSE159929/GSE159929_Bladder.RData")
load("~/Downloads/scRNA/EC/Bladder/GSE169379/GSE169379_Bladder.RData")
GSE129845_Bladder = subset(GSE129845_Bladder, DoubletFinder == "Singlet")
GSE159929_Bladder = subset(GSE159929_Bladder, DoubletFinder == "Singlet")
GSE169379_Bladder = subset(GSE169379_Bladder, DoubletFinder == "Singlet")
# integrated data
scRNA = merge(GSE129845_Bladder, y=c(GSE159929_Bladder,GSE169379_Bladder))
rm(GSE129845_Bladder,GSE159929_Bladder,GSE169379_Bladder)
save(scRNA,file="Bladder_merge.RData")
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
save(scRNA,file="Bladder_seurat.RData")
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE) 
head(Idents(scRNA), 5)
save(scRNA,file="Bladder_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","DCN","COL1A2","PECAM1","VWF","FLT1","CDH5","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)
table(scRNA$seurat_clusters)
Bladder_Endo <- scRNA[,scRNA$seurat_clusters %in% c("17","28")]
Bladder_Endo <- CreateSeuratObject(Bladder_Endo@assays$RNA@counts, meta.data = Bladder_Endo@meta.data)
save(Bladder_Endo,file="Bladder_Endo.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Brain
##########################################################################################################################
# loading data
#load data
load("~/Aging/Endo/ACC_Endo.RData")
load("~/Aging/Endo/AMY_Endo.RData")
load("~/Aging/Endo/EC_Endo.RData")
load("~/Aging/Endo/HIP_Endo.RData")
load("~/Aging/Endo/MB_Endo.RData")
load("~/Aging/Endo/PFC_Endo.RData")
load("~/Aging/Endo/Retina_Endo.RData")
load("~/Aging/Endo/SC_Endo.RData")
scRNA = merge(ACC_Endo, y=c(AMY_Endo,EC_Endo,HIP_Endo,MB_Endo,PFC_Endo,Retina_Endo,SC_Endo))
rm(ACC_Endo,AMY_Endo,EC_Endo,HIP_Endo,MB_Endo,PFC_Endo,Retina_Endo,SC_Endo)
save(scRNA,file="Brain_Endo.RData")
load("/EC/Brain/Brain_Endo.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","DCN","COL1A2","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)
table(scRNA$seurat_clusters)
Brain_Endo <- scRNA[,scRNA$Region %in% c("ACC","EC","HIP","MB","PFC","Retina","SC")]
Brain_Endo$Region <- as.factor(as.character(Brain_Endo$Region))
Brain_Endo = merge(Brain_Endo, y=c(OC_OTC_Endo))
##
current.cluster.ids <- c("ACC", "EC", "HIP", "MB", "OC", "OTC","PFC", "Retina", "SC")
new.cluster.ids <- c("Brain", "Brain", "Brain", "Brain", "Brain", "Brain","Brain", "Retina", "SC") 
Brain_Endo$Organ <- plyr::mapvalues(x = as.character(Brain_Endo$Region), from = current.cluster.ids, to = new.cluster.ids)
table(Brain_Endo$Organ,Brain_Endo$Region)
Brain_Endo <- Brain_Endo[,Brain_Endo$Region %in% c("ACC","EC","HIP","MB","OC","OTC","PFC","SC")]
Brain_Endo$Region <- as.factor(as.character(Brain_Endo$Region))
save(Brain_Endo,file="Brain_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Esophagus
##########################################################################################################################
# loading data
load("GSM4850580.RData")
GSE159929_Esophagus[["datasets"]] <- "GSE159929_Esophagus"
GSE159929_Esophagus = subset(GSE159929_Esophagus, DoubletFinder == "Singlet")
# integrated data
scRNA=GSE159929_Esophagus
DefaultAssay(scRNA) <- "SCT"
markers.to.plot <- c("MYH11","ACTA2","RGS5","DCN","COL1A2","PDGFRA","PDGFRB","PECAM1","VWF","FLT1","CDH5","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)
table(scRNA$seurat_clusters)
GSE159929_Endo <- scRNA[,scRNA$seurat_clusters %in% c("6","8")]
GSE159929_Endo <- CreateSeuratObject(GSE159929_Endo@assays$RNA@counts, meta.data = GSE159929_Endo@meta.data)
Esophagus_Endo <- GSE159929_Endo
save(Esophagus_Endo,file="Esophagus_Endo1.RData")
##
DefaultAssay(scRNA) <- "SCT"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# MG
##########################################################################################################################
# loading data
load("GSE161529_MG.RData")
GSE161529_MG = subset(GSE161529_MG, DoubletFinder == "Singlet")
scRNA=GSE161529_MG
# integrated data
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
save(scRNA,file="MG_seurat.RData")
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:30)
scRNA <- RunTSNE(scRNA, dims = 1:30,check_duplicates = FALSE)  ##耗时久
head(Idents(scRNA), 5)
save(scRNA,file="MG_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","DCN","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)
table(scRNA$seurat_clusters)
MG_Endo <- scRNA[,scRNA$seurat_clusters %in% c("7","8","13")]
MG_Endo <- CreateSeuratObject(MG_Endo@assays$RNA@counts, meta.data = MG_Endo@meta.data)
MG_Endo[["Region"]] <- "MG"
MG_Endo[["Organ"]] <- "MG"
save(MG_Endo,file="MG_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Heart
##########################################################################################################################
# loading data
load("ERP123138_Heart.RData")
load("SCP498_Heart.RData")
ERP123138_Heart = subset(ERP123138_Heart, DoubletFinder == "Singlet")
SCP498_Heart$Sample <- SCP498_Heart$SampleID
SCP498_Heart$SampleID <- paste(SCP498_Heart$Sample, SCP498_Heart$Region, sep = "_")
#SCP498_Heart = subset(SCP498_Heart, DoubletFinder == "Singlet")
# integrated data
scRNA = merge(ERP123138_Heart,y=c(SCP498_Heart))
scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
})
scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- ScaleData(x, features = scRNA.features, verbose = FALSE)
  x <- RunPCA(x, features = scRNA.features, verbose = FALSE)
})
# select reference M-96-P1681_LA;M-109-P1702_RA;F-104-P1708_LV;F-116-P1723_RV;F-5-D1_RV;F-76-H5_AX;M-40-D6_AX;M-66-H3_LV
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, reference = c(5,40,66,76,96,104,109,116), anchor.features = 3000,reduction = "rpca", dims = 1:50)
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:50)
scRNA
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)
ElbowPlot(scRNA, ndims = 50)
scRNA <- FindNeighbors(scRNA, dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:30)
scRNA <- RunTSNE(scRNA, dims = 1:30,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="scRNA_Heart_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","DCN","COL1A2","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)
table(scRNA$seurat_clusters)
Heart_Endo <- scRNA[,scRNA$seurat_clusters %in% c("8","21","23")]
Heart_Endo <- CreateSeuratObject(Heart_Endo@assays$RNA@counts, meta.data = Heart_Endo@meta.data)
save(Heart_Endo,file="Heart_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)

##########################################################################################################################
# intestinal
##########################################################################################################################
# loading data
load("E_MATB_9543_GI.RData")
scRNA = subset(E_MATB_9543_GI, DoubletFinder == "Singlet")
# integrated data
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="GI_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","ACTA2","RGS5","DCN","COL1A2","PDGFRA","PDGFRB","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)
GI_Endo <- scRNA[,scRNA$seurat_clusters %in% c("26","30")]
GI_Endo <- CreateSeuratObject(GI_Endo@assays$RNA@counts, meta.data = GI_Endo@meta.data)
current.cluster.ids <- c("appendix","ascending colon","caecum","descending colon","duodenum","ileum","jejunum","lymph node","sigmoid colon","transverse colon")
new.cluster.ids <- c("colon","colon","colon","colon","colon","small intestine","small intestine","lymph node","colon","colon") 
GI_Endo$Region <- plyr::mapvalues(x = as.character(GI_Endo$Region), from = current.cluster.ids, to = new.cluster.ids)
table(GI_Endo@meta.data$Region)
save(GI_Endo,file="GI_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Kidney
##########################################################################################################################
# loading data
load("GSE131882_Kindey.RData")
load("GSE171668_Kidney.RData")
GSE131882_Kindey = subset(GSE131882_Kindey, DoubletFinder == "Singlet")
GSE171668_Kidney = subset(GSE171668_Kidney, doublet == "FALSE")
# integrated data
scRNA = merge(GSE131882_Kindey, y=c(GSE171668_Kidney))
rm(GSE131882_Kindey,GSE171668_Kidney)
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
save(scRNA,file="Kidney_seurat.RData")
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:30)
scRNA <- RunTSNE(scRNA, dims = 1:30,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="Kidney_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","DCN","COL1A2","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)
table(scRNA$seurat_clusters)
Kidney_Endo <- scRNA[,scRNA$seurat_clusters %in% c("8","21","23")]
Kidney_Endo <- CreateSeuratObject(Kidney_Endo@assays$RNA@counts, meta.data = Kidney_Endo@meta.data)
Kidney_Endo[["Region"]] <- "Kidney"
Kidney_Endo[["Organ"]] <- "Kidney"
save(Kidney_Endo,file="Kidney_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Liver
##########################################################################################################################
# loading data
load("GSE115469_Liver.RData")
load("GSE136103_Liver.RData")
load("GSE158723_Liver.RData")
load("GSE185477_Liver.RData")
GSE115469_Liver = subset(GSE115469_Liver, DoubletFinder == "Singlet")
GSE136103_Liver = subset(GSE136103_Liver, DoubletFinder == "Singlet")
GSE158723_Liver = subset(GSE158723_Liver, DoubletFinder == "Singlet")
GSE185477_Liver = subset(GSE185477_Liver, DoubletFinder == "Singlet")
# integrated data
scRNA = merge(GSE115469_Liver, y=c(GSE136103_Liver,GSE158723_Liver,GSE185477_Liver))
rm(GSE115469_Liver,GSE136103_Liver,GSE158723_Liver,GSE185477_Liver)
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="Liver_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","ACTA2","RGS5","DCN","COL1A2","PDGFRA","PDGFRB","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)
table(scRNA$seurat_clusters)
Liver_Endo <- scRNA[,scRNA$seurat_clusters %in% c("4","16","23","24","28")]
Liver_Endo <- CreateSeuratObject(Liver_Endo@assays$RNA@counts, meta.data = Liver_Endo@meta.data)
Liver_Endo = merge(Liver_Endo, y=c(E_MATB_10553_Endo))
Liver_Endo[["Region"]] <- "Liver"
Liver_Endo[["Organ"]] <- "Liver"
save(Liver_Endo,file="Liver_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Liver
##########################################################################################################################
# loading data
load("GSE135893_Lung.RData")
load("GSE136831_Lung.RData")
load("GSE164829_Lung.RData")
load("PRJEB31843_Lung.RData")
GSE135893_Lung = subset(GSE135893_Lung, DoubletFinder == "Singlet")
GSE136831_Lung = subset(GSE136831_Lung, DoubletFinder == "Singlet")
# integrated data
scRNA = merge(GSE135893_Lung, y=c(GSE136831_Lung,GSE164829_Lung,PRJEB31843_Lung))
rm(GSE135893_Lung,GSE136831_Lung,GSE164829_Lung,PRJEB31843_Lung)
save(scRNA,file="Lung_merge.RData")
scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
})
scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- ScaleData(x, features = scRNA.features, verbose = FALSE)
  x <- RunPCA(x, features = scRNA.features, verbose = FALSE)
})
# select reference
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, reference = c(1,4,11,15,41,44), anchor.features = 3000,reduction = "rpca", dims = 1:50)
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="Lung_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","ACTA2","RGS5","DCN","COL1A2","PDGFRA","PDGFRB","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)
table(scRNA$seurat_clusters)
Skin_Endo <- scRNA[,scRNA$seurat_clusters %in% c("11","21","23")]
Skin_Endo <- CreateSeuratObject(Skin_Endo@assays$RNA@counts, meta.data = Skin_Endo@meta.data)
save(Skin_Endo,file="Skin_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Decidua
##########################################################################################################################
# loading data
load("E_MTAB_6701_Decidua.RData")
E_MTAB_6701_Decidua = subset(E_MTAB_6701_Decidua, DoubletFinder == "Singlet")
scRNA = E_MTAB_6701_Decidua
# integrated data
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="Decidua_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","ACTA2","RGS5","DCN","COL1A2","PDGFRA","PDGFRB","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)
table(scRNA$seurat_clusters)
Decidua_Endo <- scRNA[,scRNA$seurat_clusters %in% c("9","14","20")]
Decidua_Endo <- CreateSeuratObject(Decidua_Endo@assays$RNA@counts, meta.data = Decidua_Endo@meta.data)
save(Decidua_Endo,file="Decidua_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Placenta
##########################################################################################################################
# loading data
load("E_MTAB_6701_Placenta.RData")
E_MTAB_6701_Placenta = subset(E_MTAB_6701_Placenta, DoubletFinder == "Singlet")
scRNA = E_MTAB_6701_Placenta
# integrated data
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
save(scRNA,file="Placenta_seurat.RData")
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="Placenta_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","ACTA2","RGS5","DCN","COL1A2","PDGFRA","PDGFRB","PECAM1","VWF","FLT1","CDH5","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
markers.to.plot <- c("DLK1","EGFR","GATA3","PAGE4","LGALS13","HLA-G")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)
table(scRNA$seurat_clusters)
Placenta_Endo <- scRNA[,scRNA$seurat_clusters %in% c("23")]
Placenta_Endo <- CreateSeuratObject(Placenta_Endo@assays$RNA@counts, meta.data = Placenta_Endo@meta.data)
save(Placenta_Endo,file="Placenta_Endo.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Muscle
##########################################################################################################################
# loading data
load("GSE143704_Muscle.RData")
load("GSE214544_Muscle.RData")
load("GSE159929_Muscle.RData")
GSE143704_Muscle = subset(GSE143704_Muscle, DoubletFinder == "Singlet")
GSE214544_Muscle = subset(GSE214544_Muscle, DoubletFinder == "Singlet")
GSE159929_Muscle = subset(GSE159929_Muscle, DoubletFinder == "Singlet")
# integrated data
scRNA = merge(GSE143704_Muscle, y=c(GSE214544_Muscle,GSE159929_Muscle))
rm(GSE143704_Muscle,GSE214544_Muscle,GSE159929_Muscle)
save(scRNA,file="Muscle_merge.RData")
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
save(scRNA,file="Muscle_seurat.RData")
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="Muscle_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","DCN","COL1A2","PECAM1","VWF","FLT1","CDH5","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)
table(scRNA$seurat_clusters)
Muscle_Endo <- scRNA[,scRNA$seurat_clusters %in% c("0","3","7","9","10","12","14","19","23","31")]
Muscle_Endo <- CreateSeuratObject(Muscle_Endo@assays$RNA@counts, meta.data = Muscle_Endo@meta.data)
save(Muscle_Endo,file="Muscle_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Ovary
##########################################################################################################################
# loading data
load("GSE118127_Ovary.RData")
load("GSE213216_Ovary.RData")
GSE118127_Ovary = subset(GSE118127_Ovary, DoubletFinder == "Singlet")
GSE213216_Ovary = subset(GSE213216_Ovary, DoubletFinder == "Singlet")
# integrated data
scRNA = merge(GSE118127_Ovary, y=c(GSE213216_Ovary))
rm(GSE118127_Ovary,GSE213216_Ovary)
save(scRNA,file="Ovary_merge.RData")
scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
})
scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- ScaleData(x, features = scRNA.features, verbose = FALSE)
  x <- RunPCA(x, features = scRNA.features, verbose = FALSE)
})

# select reference
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, reference = c(1,2,3,4,30,31,32,33), anchor.features = 3000,reduction = "rpca", dims = 1:50)
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE) 
head(Idents(scRNA), 5)
save(scRNA,file="Ovary_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","ACTA2","RGS5","DCN","COL1A2","PDGFRA","PDGFRB","PECAM1","VWF","FLT1","CLDN5","CDH5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)
table(scRNA$seurat_clusters)
Ovary_Endo <- scRNA[,scRNA$seurat_clusters %in% c("7","10","12","18")]
Ovary_Endo <- CreateSeuratObject(Ovary_Endo@assays$RNA@counts, meta.data = Ovary_Endo@meta.data)
save(Ovary_Endo,file="Ovary_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Pancreas
##########################################################################################################################
# identified EC
adult_pancreas_2020 <- readRDS("adult_pancreas_2020.rds")
scRNA <- adult_pancreas_2020
scRNA$SampleID <- scRNA$sample_ID
scRNA$Age <- scRNA$age
scRNA$Sex <- scRNA$sex
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","DCN","COL1A2","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
table(scRNA$Cluster)
Pancreas_Endo <- scRNA[,scRNA$Cluster %in% c("Endothelial")]
Pancreas_Endo <- Pancreas_Endo[,Pancreas_Endo$patient_ID %in% c("AFES365","AFES448","TUM_13","TUM_25","TUM_C1")]
Pancreas_Endo <- CreateSeuratObject(Pancreas_Endo@assays$RNA@counts, meta.data = Pancreas_Endo@meta.data)
Pancreas_Endo[["Region"]] <- "Pancreas"
Pancreas_Endo[["Organ"]] <- "Pancreas"
save(Pancreas_Endo,file="Pancreas_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# PIT
##########################################################################################################################
# loading data
load("/GSE178454_PIT/GSE178454_PIT.RData")
GSE178454_PIT = subset(GSE178454_PIT, DoubletFinder == "Singlet")
# integrated data
scRNA = GSE178454_PIT
rm(GSE178454_PIT)
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:50)
scRNA
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 100, verbose = FALSE)
ElbowPlot(scRNA, ndims = 50)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 1.2)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="PIT_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","DCN","COL1A2","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)
table(scRNA$seurat_clusters)
PIT_Endo <- scRNA[,scRNA$seurat_clusters %in% c("23","25")]
PIT_Endo <- CreateSeuratObject(PIT_Endo@assays$RNA@counts, meta.data = PIT_Endo@meta.data)
PIT_Endo[["Organ"]] <- "PIT"
save(PIT_Endo,file="PIT_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Skin
##########################################################################################################################
# loading data
load("HRA000395_Skin.RData")
load("PRJCA002692_Skin.RData")
load("GSE159929_Skin.RData")
HRA000395_Skin = subset(HRA000395_Skin, DoubletFinder == "Singlet")
PRJCA002692_Skin = subset(PRJCA002692_Skin, DoubletFinder == "Singlet")
GSE159929_Skin = subset(GSE159929_Skin, DoubletFinder == "Singlet")
# integrated data
scRNA = merge(HRA000395_Skin, y=c(PRJCA002692_Skin,GSE159929_Skin))
rm(HRA000395_Skin,PRJCA002692_Skin,GSE159929_Skin)
save(scRNA,file="skin_merge.RData")
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}

# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
save(scRNA,file="skin_seurat.RData")
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="skin_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","ACTA2","RGS5","DCN","COL1A2","PDGFRA","PDGFRB","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)
table(scRNA$seurat_clusters)
Skin_Endo <- scRNA[,scRNA$seurat_clusters %in% c("0","1","6","7","8","11","13","15","18","26")]
Skin_Endo <- CreateSeuratObject(Skin_Endo@assays$RNA@counts, meta.data = Skin_Endo@meta.data)
save(Skin_Endo,file="Skin_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Testis
##########################################################################################################################
# loading data
load("GSE149512_Testis.RData")
load("GSE182786_Testis.RData")
GSE149512_Testis = subset(GSE149512_Testis, DoubletFinder == "Singlet")
GSE182786_Testis = subset(GSE182786_Testis, DoubletFinder == "Singlet")
# integrated data
scRNA = merge(GSE149512_Testis, y=c(GSE182786_Testis))
rm(GSE149512_Testis,GSE182786_Testis)
save(scRNA,file="Testis_merge.RData")
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
save(scRNA,file="Testis_seurat.RData")
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE) 
head(Idents(scRNA), 5)
save(scRNA,file="Testis_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","ACTA2","RGS5","DCN","COL1A2","PDGFRA","PDGFRB","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)
table(scRNA$seurat_clusters)
Testis_Endo <- scRNA[,scRNA$seurat_clusters %in% c("11","18","20","31")]
Testis_Endo <- CreateSeuratObject(Testis_Endo@assays$RNA@counts, meta.data = Testis_Endo@meta.data)
save(Testis_Endo,file="Testis_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Thyroid
##########################################################################################################################
# loading data
load("GSE182416_Thyroid.RData")
load("GSE191288_Thyroid.RData")
GSE191288_Thyroid = subset(GSE191288_Thyroid, DoubletFinder == "Singlet")
# integrated data
scRNA = merge(GSE182416_Thyroid, y=c(GSE191288_Thyroid))
rm(GSE182416_Thyroid,GSE191288_Thyroid)
save(scRNA,file="Thyroid_merge.RData")
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
save(scRNA,file="Thyroid_seurat.RData")
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:30)
scRNA <- RunTSNE(scRNA, dims = 1:30,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="Thyroid_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","DCN","COL1A2","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
table(scRNA$seurat_clusters)
Thyroid_Endo <- scRNA[,scRNA$seurat_clusters %in% c("8","17","18","21","23")]
Thyroid_Endo <- CreateSeuratObject(Thyroid_Endo@assays$RNA@counts, meta.data = Thyroid_Endo@meta.data)
Thyroid_Endo[["Organ"]] <- "Thyroid"
save(Thyroid_Endo,file="Thyroid_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Trachea
##########################################################################################################################
# loading data
load("GSE159929_Trachea.RData")
# identified EC
GSE159929_Trachea = subset(GSE159929_Trachea, DoubletFinder == "Singlet")
scRNA=GSE159929_Trachea
DefaultAssay(scRNA) <- "SCT"
markers.to.plot <- c("MYH11","ACTA2","RGS5","DCN","COL1A2","PDGFRA","PDGFRB","PECAM1","VWF","FLT1","CDH5","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)
table(scRNA$seurat_clusters)
Trachea_Endo <- scRNA[,scRNA$seurat_clusters %in% c("2")]
Trachea_Endo <- CreateSeuratObject(Trachea_Endo@assays$RNA@counts, meta.data = Trachea_Endo@meta.data)
save(Trachea_Endo,file="Trachea_Endo1.RData")
##
DefaultAssay(scRNA) <- "SCT"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
##########################################################################################################################
# Uterus
##########################################################################################################################
# loading data
load("E_MTAB_10287_Uterus.RData")
load("GSE213216_Uterus.RData")
E_MTAB_10287_Uterus = subset(E_MTAB_10287_Uterus, DoubletFinder == "Singlet")
GSE213216_Uterus = subset(GSE213216_Uterus, DoubletFinder == "Singlet")
# integrated data
scRNA = merge(E_MTAB_10287_Uterus, y=c(GSE213216_Uterus))
rm(E_MTAB_10287_Uterus,GSE213216_Uterus)
save(scRNA,file="Uterus_merge.RData")
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
# this command creates an 'integrated' data assay
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
save(scRNA,file="Uterus_seurat.RData")
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 30, verbose = FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="Uterus_seurat.RData")
# identified EC
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MYH11","DCN","COL1A2","PECAM1","VWF","FLT1","CLDN5","IFI27","AQP1","SPRY1","MT1M","PTPRC","LYZ","CD3D","KLRD1","CD79A","MZB1","CPA3")
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
table(scRNA$seurat_clusters)
Thyroid_Endo <- scRNA[,scRNA$seurat_clusters %in% c("8","17","18","21","23")]
Thyroid_Endo <- CreateSeuratObject(Thyroid_Endo@assays$RNA@counts, meta.data = Thyroid_Endo@meta.data)
Thyroid_Endo[["Organ"]] <- "Thyroid"
save(Thyroid_Endo,file="Thyroid_Endo1.RData")
##
DefaultAssay(scRNA) <- "RNA"
celltype_umap<- FeaturePlot(scRNA, features = c("PECAM1","VWF"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo1_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("FLT1","CLDN5"),ncol = 1,pt.size = 0.001,cols =c("lightgrey","#923331"),raster=FALSE) 
ggsave(filename = "Endo2_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)



