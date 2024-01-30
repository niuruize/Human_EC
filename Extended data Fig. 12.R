
################################################################################
# Extended data Fig. 12
################################################################################
# The difference in cell numbers between female and male for ECs vascular bed subtypes in different tissues was calculated separately
#devtools::install_github("MarioniLab/miloR", red = "devel")
library(miloR);library(Seurat);library(ggplot2);library(SingleCellExperiment);library(SeuratWrappers);library(scuttle);
library(scater);library(scales);library(forcats);library(data.table);library(stringr);library(dplyr);library(scran);library(patchwork)
unique(scRNA$Sex)
table(scRNA$Sex)
load('Brain_Harmony.RData')
scRNA = scRNA[,scRNA$Sex %in% c("F","M")]
scRNA$Sex <- as.factor(as.character(scRNA$Sex))
DefaultAssay(scRNA)="RNA"
scRNA@assays$RNA@data <- scRNA@assays$RNA@counts
scRNA <- NormalizeData(scRNA)
#--Define cell neighbourhoods
scRNA1 <- as.SingleCellExperiment(scRNA)
traj_milo <- miloR::Milo(scRNA1)
reducedDim(traj_milo, "UMAP") <- scRNA@reductions$umap@cell.embeddings
reducedDim(traj_milo, "PCA") <- scRNA@reductions$pca@cell.embeddings
reducedDim(traj_milo, "TSNE") <- scRNA@reductions$tsne@cell.embeddings

plotUMAP(traj_milo, colour_by="celltype") + plotUMAP(traj_milo, colour_by="Sex")
traj_milo <- buildGraph(traj_milo, k = 10, d = 10)
traj_milo <- makeNhoods(traj_milo, prop = 0.2, k = 10, d=30, refined = TRUE, refinement_scheme="graph")
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), samples="SampleID")
#--Differential abundance testing with contrasts
traj_design <- data.frame(colData(traj_milo))[,c("SampleID", "Sex")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$SampleID
## Reorder rownames to match columns of nhoodCounts(milo)
traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]
table(traj_design$Sex)
rownames(traj_design) <- traj_design$SampleID
contrast.1 <- c("SexF - SexM")
da_results <- testNhoods(traj_milo, design = ~ 0 + Sex, design.df = traj_design, model.contrasts = contrast.1, fdr.weighting="graph-overlap")
da_results %>% arrange(- SpatialFDR) %>% head() 
table(da_results$SpatialFDR < 0.1)

#Visualize neighbourhoods displaying DA
traj_milo <- buildNhoodGraph(traj_milo)
plotUMAP(traj_milo, colour_by = "celltype") + plotNhoodGraphDA(traj_milo, da_results, alpha=0.05) + plot_layout(guides="collect")

#可视化1----------------------------------------------------------------------------
#UMAP图展示neighborhoods
pdf('milo_PMN_Sex.pdf', width=18, height=8)
plotUMAP(traj_milo, colour_by = "celltype",point_size = 2) + theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
  plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20)) +
  plotNhoodGraphDA(traj_milo, da_results, alpha=0.1) +
  scale_fill_gradient2(low="#070091",#修改颜色
                       mid="lightgrey",
                       high="#910000", 
                       name="log2FC",
                       limits=c(-5,5),
                       oob=squish) + 
  plot_layout(guides="collect")
dev.off()
#可视化2----------------------------------------------------------------------------
#蜂群图展示celltype logFC变化
da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "celltype")
pdf('milo_beeswarm_Sex.pdf', width=3, height=4)
plotDAbeeswarm(da_results, group.by = "celltype") +
  scale_color_gradient2(low="#070091",
                        mid="lightgrey",
                        high="#910000",
                        #limits=c(-5,5),
                        oob=squish) +
  labs(x="", y="Log2 Fold Change") +
  theme_bw(base_size=10)+
  theme(axis.text = element_text(colour = 'black')) 
dev.off()





