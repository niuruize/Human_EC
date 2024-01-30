
################################################################################
# Extended data Fig. 9A
################################################################################
# The difference in cell numbers between young and old for ECs vascular bed subtypes in different tissues was calculated separately
#devtools::install_github("MarioniLab/miloR", red = "devel")
library(miloR);library(Seurat);library(ggplot2);library(SingleCellExperiment);library(SeuratWrappers);library(scuttle);
library(scater);library(scales);library(forcats);library(data.table);library(stringr);library(dplyr);library(scran);library(patchwork)
unique(scRNA$Group)
table(scRNA$Group)
scRNA = scRNA[,scRNA$Group %in% c("Adult","Aging")]
scRNA$Group <- as.factor(as.character(scRNA$Group))
DefaultAssay(scRNA)="RNA"
scRNA@assays$RNA@data <- scRNA@assays$RNA@counts
scRNA <- NormalizeData(scRNA)
scRNA1 <- as.SingleCellExperiment(scRNA)
#--Define cell neighbourhoods
scRNA1 <- as.SingleCellExperiment(scRNA)
traj_milo <- miloR::Milo(scRNA1)
reducedDim(traj_milo, "UMAP") <- scRNA@reductions$umap@cell.embeddings
reducedDim(traj_milo, "PCA") <- scRNA@reductions$pca@cell.embeddings
reducedDim(traj_milo, "TSNE") <- scRNA@reductions$tsne@cell.embeddings

plotUMAP(traj_milo, colour_by="celltype") + plotUMAP(traj_milo, colour_by="Group")
traj_milo <- buildGraph(traj_milo, k = 10, d = 10)
traj_milo <- makeNhoods(traj_milo, prop = 0.2, k = 10, d=30, refined = TRUE, refinement_scheme="graph")
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), samples="SampleID")
#--Differential abundance testing with contrasts
traj_design <- data.frame(colData(traj_milo))[,c("SampleID", "Group")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$SampleID
## Reorder rownames to match columns of nhoodCounts(milo)
traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]
table(traj_design$Group)
rownames(traj_design) <- traj_design$SampleID
contrast.1 <- c("GroupAging - GroupAdult")
da_results <- testNhoods(traj_milo, design = ~ 0 + Group, design.df = traj_design, model.contrasts = contrast.1, fdr.weighting="graph-overlap")
da_results %>% arrange(- SpatialFDR) %>% head() 
table(da_results$SpatialFDR < 0.1)

#Visualize neighbourhoods displaying DA
traj_milo <- buildNhoodGraph(traj_milo)
plotUMAP(traj_milo, colour_by = "celltype") + plotNhoodGraphDA(traj_milo, da_results, alpha=0.05) + plot_layout(guides="collect")

#
#UMAP
pdf('milo_PMN_Group.pdf', width=18, height=8)
plotUMAP(traj_milo, colour_by = "celltype",point_size = 2) + theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
  plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20)) +
  plotNhoodGraphDA(traj_milo, da_results, alpha=0.1) +
  scale_fill_gradient2(low="#070091",
                       mid="lightgrey",
                       high="#910000", 
                       name="log2FC",
                       limits=c(-5,5),
                       oob=squish) + 
  plot_layout(guides="collect")
dev.off()
#
da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "celltype")
pdf('milo_beeswarm_Group.pdf', width=3, height=4)
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

################################################################################
# Extended data Fig. 9B
################################################################################
# Aorta
library(Seurat);library(tidyverse);library(cowplot);library(Matrix);library(readxl);library(ggpubr)
load('Aorta_EC.RData')
Idents(scRNA)="Group"
table(scRNA$Group)
Interferon_gene<-read_xlsx("aging.xlsx")
gene<-as.list(Interferon_gene)
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "aging")
colnames(AB@meta.data)[29]<-"Interferon_Score"
##
b<-FetchData(AB, vars = c("Group","Interferon_Score"))
b[["Group"]]<-factor(b[["Group"]], levels=c("Adult","Aging"))
my_comparisons <- list(c("Adult", "Aging"))
ggviolin(b, x = "Group", y = "Interferon_Score", fill = "Group",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("Interferon_7.pdf",width = 10,height = 15,units = "cm")

# Liver
library(Seurat);library(tidyverse);library(cowplot);library(Matrix);library(readxl);library(ggpubr)
load('Liver_EC.RData')
Idents(scRNA)="Group"
table(scRNA$Group)
Interferon_gene<-read_xlsx("aging.xlsx")
gene<-as.list(Interferon_gene)
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "aging")
colnames(AB@meta.data)[28]<-"Interferon_Score"
##
b<-FetchData(AB, vars = c("Group","Interferon_Score"))
b[["Group"]]<-factor(b[["Group"]], levels=c("Adult","Aging"))
my_comparisons <- list(c("Adult", "Aging"))
ggboxplot(b, x = "Group", y = "Interferon_Score",combine = TRUE,add = "jitter", 
          add.params = list(size=0.1, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=9, face="italic"), repel = TRUE,
          color = "Group", palette = "npg")+stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=45,size=6))
ggsave("Interferon_6.pdf",width = 10,height = 15,units = "cm")

ggviolin(b, x = "Group", y = "Interferon_Score", fill = "Group",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("Interferon_7.pdf",width = 10,height = 15,units = "cm")

################################################################################
# Extended data Fig. 9C
################################################################################
# heart
library(Seurat);library(tidyverse);library(cowplot);library(Matrix);library(readxl);library(ggpubr)
load('GSE183852_heart_EC.RData')
Human_aging_up_gene<-read_xlsx("aging.xlsx")
gene<-as.list(Human_aging_up_gene)
Idents(scRNA)="Names"
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "aging")
colnames(AB@meta.data)[57]<-"Human_aging_up_Score"
#
P1<-VlnPlot(AB, features = 'Human_aging_up_Score')
ggsave(filename = "Human_aging_up_1.pdf", plot = P1, device = 'pdf', width = 20, height = 15, units = 'cm')
rm('P1')  
#
P1<-RidgePlot(AB, features = 'Human_aging_up_Score', ncol = 1) 
ggsave(filename = "Human_aging_up_2.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1')  
#
data<-FetchData(AB, vars = c("Names","Human_aging_up_Score"))
P1<-ggplot(data,aes(Names, Human_aging_up_Score))+
  geom_boxplot()+theme_bw()+RotatedAxis()
ggsave(filename = "Human_aging_up_3.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 
##Distribution
b<-FetchData(AB, vars = c("Group","Human_aging_up_Score"))
b[["Group"]]<-factor(b[["Group"]], levels=c("Young","Old"))
P1<-ggdensity(b, x = "Human_aging_up_Score",
          add = "mean", rug = TRUE,
          color = "Group", fill = "Group",
          palette = c("#00AFBB", "#E7B800"))
ggsave(filename = "Human_aging_up_4.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 
P1<-gghistogram(b, x = "Human_aging_up_Score",bins = 30, add = "mean", rug = TRUE,
            color = "Group", fill = "Group",
            palette = c("#00AFBB", "#E7B800"))
ggsave(filename = "Human_aging_up_5.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 
##
b<-FetchData(AB, vars = c("Group","Human_aging_up_Score"))
b[["Group"]]<-factor(b[["Group"]], levels=c("Young","Middle","Old"))
my_comparisons <- list(c("Young", "Old"))
ggboxplot(b, x = "Group", y = "Human_aging_up_Score",combine = TRUE,add = "jitter", 
          add.params = list(size=0.1, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=9, face="italic"), repel = TRUE,
          color = "Group", palette = "npg")+stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=45,size=6))
ggsave("Human_aging_up_6.pdf",width = 10,height = 15,units = "cm")

ggviolin(b, x = "Group", y = "Human_aging_up_Score", fill = "Group",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("Human_aging_up_7.pdf",width = 10,height = 15,units = "cm")

# lung
load('GSE171668_Lung_EC.RData')
Human_aging_up_gene<-read_xlsx("aging.xlsx")
gene<-as.list(Human_aging_up_gene)
Idents(scRNA)="Names"
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "aging")
colnames(AB@meta.data)[31]<-"Human_aging_up_Score"
#
P1<-VlnPlot(AB, features = 'Human_aging_up_Score')
ggsave(filename = "Human_aging_up_1.pdf", plot = P1, device = 'pdf', width = 20, height = 15, units = 'cm')
rm('P1')  
#
P1<-RidgePlot(AB, features = 'Human_aging_up_Score', ncol = 1) 
ggsave(filename = "Human_aging_up_2.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1')  
#
data<-FetchData(AB, vars = c("Names","Human_aging_up_Score"))
P1<-ggplot(data,aes(Names, Human_aging_up_Score))+
  geom_boxplot()+theme_bw()+RotatedAxis()
ggsave(filename = "Human_aging_up_3.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 
##Distribution
b<-FetchData(AB, vars = c("Group","Human_aging_up_Score"))
b[["Group"]]<-factor(b[["Group"]], levels=c("Young","Old"))
P1<-ggdensity(b, x = "Human_aging_up_Score",
          add = "mean", rug = TRUE,
          color = "Group", fill = "Group",
          palette = c("#00AFBB", "#E7B800"))
ggsave(filename = "Human_aging_up_4.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 
P1<-gghistogram(b, x = "Human_aging_up_Score",bins = 30, add = "mean", rug = TRUE,
            color = "Group", fill = "Group",
            palette = c("#00AFBB", "#E7B800"))
ggsave(filename = "Human_aging_up_5.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 
##
b<-FetchData(AB, vars = c("Group","Human_aging_up_Score"))
b[["Group"]]<-factor(b[["Group"]], levels=c("Young","Middle","Old"))
my_comparisons <- list(c("Young", "Old"))
ggboxplot(b, x = "Group", y = "Human_aging_up_Score",combine = TRUE,add = "jitter", 
          add.params = list(size=0.1, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=9, face="italic"), repel = TRUE,
          color = "Group", palette = "npg")+stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=45,size=6))
ggsave("Human_aging_up_6.pdf",width = 10,height = 15,units = "cm")

ggviolin(b, x = "Group", y = "Human_aging_up_Score", fill = "Group",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("Human_aging_up_7.pdf",width = 10,height = 15,units = "cm")

# liver
load('GSE171668_Liver_EC.RData')
Human_aging_up_gene<-read_xlsx("aging.xlsx")
gene<-as.list(Human_aging_up_gene)
Idents(scRNA)="Names"
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "aging")
colnames(AB@meta.data)[27]<-"Human_aging_up_Score"
#
P1<-VlnPlot(AB, features = 'Human_aging_up_Score')
ggsave(filename = "Human_aging_up_1.pdf", plot = P1, device = 'pdf', width = 20, height = 15, units = 'cm')
rm('P1')  
#
P1<-RidgePlot(AB, features = 'Human_aging_up_Score', ncol = 1) 
ggsave(filename = "Human_aging_up_2.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1')  
#
data<-FetchData(AB, vars = c("Names","Human_aging_up_Score"))
P1<-ggplot(data,aes(Names, Human_aging_up_Score))+
  geom_boxplot()+theme_bw()+RotatedAxis()
ggsave(filename = "Human_aging_up_3.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 
##Distribution
b<-FetchData(AB, vars = c("Group","Human_aging_up_Score"))
b[["Group"]]<-factor(b[["Group"]], levels=c("Young","Old"))
P1<-ggdensity(b, x = "Human_aging_up_Score",
          add = "mean", rug = TRUE,
          color = "Group", fill = "Group",
          palette = c("#00AFBB", "#E7B800"))
ggsave(filename = "Human_aging_up_4.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 
P1<-gghistogram(b, x = "Human_aging_up_Score",bins = 30, add = "mean", rug = TRUE,
            color = "Group", fill = "Group",
            palette = c("#00AFBB", "#E7B800"))
ggsave(filename = "Human_aging_up_5.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 
##
b<-FetchData(AB, vars = c("Group","Human_aging_up_Score"))
b[["Group"]]<-factor(b[["Group"]], levels=c("Young","Middle","Old"))
my_comparisons <- list(c("Young", "Old"))
ggboxplot(b, x = "Group", y = "Human_aging_up_Score",combine = TRUE,add = "jitter", 
          add.params = list(size=0.1, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=9, face="italic"), repel = TRUE,
          color = "Group", palette = "npg")+stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=45,size=6))
ggsave("Human_aging_up_6.pdf",width = 10,height = 15,units = "cm")

ggviolin(b, x = "Group", y = "Human_aging_up_Score", fill = "Group",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("Human_aging_up_7.pdf",width = 10,height = 15,units = "cm")




