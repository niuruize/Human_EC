
################################################################################
# Extended data Fig.13A
################################################################################
library(Seurat);library(ggplot2);library(scRNAtoolVis)
load("Human_EC.RData")
p1 = jjDotPlot(scRNA, gene = c("IFITM1","IFITM2","IFITM3","IFIT1","IFIT2","IFIT3"),id = 'Organ',xtree = T,ytree = T,rescale = T,rescale.min = 0,rescale.max = 1)
ggsave(filename = "EC_Human_1.pdf", plot = p1, device = 'pdf', width = 24, height = 24, units = 'cm')
load("MacFas_EC.RData")
p1 = jjDotPlot(scRNA, gene = c("ITIM1","ITIM2","IFIT1","IFIT2","IFIT3"),id = 'Organ',xtree = T,ytree = T,rescale = T,rescale.min = 0,rescale.max = 1)
ggsave(filename = "EC_macaque_1.pdf", plot = p1, device = 'pdf', width = 24, height = 24, units = 'cm')
load("Mouse_EC.RData")
p1 = jjDotPlot(scRNA, gene = c("Ifitm1","Ifitm2","Ifitm3","Ifit1","Ifit2","Ifit3"),id = 'Organ',xtree = T,ytree = T,rescale = T,rescale.min = 0,rescale.max = 1)
ggsave(filename = "EC_Human_1.pdf", plot = p1, device = 'pdf', width = 24, height = 24, units = 'cm')
################################################################################
# Extended data Fig.13B
################################################################################
# human
library(Seurat);library(tidyverse);library(cowplot);library(Matrix);library(readxl);library(ggpubr)
load("Human_EC.RData")
Interferon_gene<-read_xlsx("aging.xlsx")
gene<-as.list(Interferon_gene)
Idents(scRNA)="Organ"
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "aging")
colnames(AB@meta.data)[17]<-"Interferon_Score"
##
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
#
b<-FetchData(AB, vars = c("Organ","Interferon_Score"))
ggboxplot(b, x = "Organ", y = "Interferon_Score",combine = TRUE,add = "jitter", 
          add.params = list(size=0.01, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=16, face="italic"), repel = TRUE,
          color = "Organ", palette = my36colors)+
  theme(axis.text.x=element_text(vjust = 1, hjust = 1, angle=90,size=8))
ggsave("Interferon_20.pdf",width = 8,height = 10,units = "cm")
# macaque
library(Seurat);library(tidyverse);library(cowplot);library(Matrix);library(readxl);library(ggpubr)
load("MacFas_EC.RData")
Interferon_gene<-read_xlsx("aging.xlsx")
gene<-as.list(Interferon_gene)
Idents(scRNA)="Organ"
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "aging")
colnames(AB@meta.data)[17]<-"Interferon_Score"
##
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
#
b<-FetchData(AB, vars = c("Organ","Interferon_Score"))
ggboxplot(b, x = "Organ", y = "Interferon_Score",combine = TRUE,add = "jitter", 
          add.params = list(size=0.01, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=16, face="italic"), repel = TRUE,
          color = "Organ", palette = my36colors)+
  theme(axis.text.x=element_text(vjust = 1, hjust = 1, angle=90,size=8))
ggsave("Interferon_20.pdf",width = 8,height = 10,units = "cm")

# mouse
load("Mouse_EC.RData")
library(Seurat);library(tidyverse);library(cowplot);library(Matrix);library(readxl);library(ggpubr)
Interferon_gene<-read_xlsx("aging.xlsx")
gene <- homologene(Interferon_gene$gene, inTax = 9606, outTax = 10090)
write.csv(gene,"gene.csv") #保存结果
Interferon_gene <- read_xlsx("aging.xlsx")
gene <- as.list(Interferon_gene)
Idents(scRNA)="Organ"
table(scRNA$Organ)
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "aging")
colnames(AB@meta.data)[12]<-"Interferon_Score"
##
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
#
b<-FetchData(AB, vars = c("Organ","Interferon_Score"))
ggboxplot(b, x = "Organ", y = "Interferon_Score",combine = TRUE,add = "jitter", 
          add.params = list(size=0.01, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=16, face="italic"), repel = TRUE,
          color = "Organ", palette = my36colors)+
  theme(axis.text.x=element_text(vjust = 1, hjust = 1, angle=90,size=8))
ggsave("Interferon_20.pdf",width = 8,height = 10,units = "cm")



################################################################################
# Extended data Fig.13C
################################################################################
load("Hu_Chi_Mac.RData")
scRNA$Species_celltype <- paste(scRNA$Species, scRNA$celltype, sep = "_")
p1 = jjDotPlot(scRNA, gene = c("IFITM1","IFITM2","IFITM3","IFIT1","IFIT2","IFIT3"),id = 'Species_celltype',xtree = T,ytree = T,rescale = T,rescale.min = 0,rescale.max = 1)
ggsave(filename = "Interferon_marker_heatmap_1.pdf", plot = p1, device = 'pdf', width = 15, height = 24, units = 'cm')



