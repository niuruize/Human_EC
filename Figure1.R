#loading R packages
library(Seurat);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

################################################################################
###---Figure 1B
################################################################################
# The input data were obtained from cross "EC_integration.R"
load("Human_EC.RData")
table(scRNA$Organ)
# plot the bar graph using Prism 9

################################################################################
###---Figure 1C
################################################################################
load("Human_EC.RData")
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175','#985a38', '#aac6d2')
p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "Organ", pt.size=0.01, label = F,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "tsne_Organ_1.pdf", plot = p1, device = 'pdf', width = 22, height = 16, units = 'cm')
ggsave(filename = "tsne_Organ_1.png", plot = p1, device = 'png', width = 22, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "Region", pt.size=0.01, label = F,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "tsne_Region_1.pdf", plot = p1, device = 'pdf', width = 24, height = 16, units = 'cm')
ggsave(filename = "tsne_Region_1.png", plot = p1, device = 'png', width = 24, height = 16, units = 'cm')

################################################################################
###---Figure 1D
################################################################################
# The input data were obtained from cross "species integration.R"
cluster.averages <- AverageExpression(scRNA,group.by = "Organ")
B.exp <- cluster.averages[["SCT"]]
hc <- as.dendrogram(hclust(dist(t(B.exp)), "ave"))
hc <- hc %>% color_branches(k = 2) %>% color_labels(k = 38)
# 
pdf("cluster_Organ.pdf", width = 6,height = 3)
plot(hc)
dev.off()
# 
pdf("circ_cluster_Organ.pdf", width = 6,height = 6)
circlize_dendrogram(hc,labels_track_height = NA,dend_track_height = 0.3)
dev.off()
################################################################################
###---Figure 1E
################################################################################
markers <- read.table("markers.txt",header = T)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
# tartget gene
annoGene <- read_xlsx("annogene.xlsx")
pdf("heatmap_marker_1.pdf", width = 6,height = 6)
AverageHeatmap(scRNA, markerGene = top10$gene,clusterAnnoName = F, showRowNames = F, markGenes = annoGene$gene)
dev.off()












