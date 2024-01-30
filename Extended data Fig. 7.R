
################################################################################
# Extended data Fig. 7A
################################################################################
library(Seurat)
# GSE188545_MTG
DefaultAssay(scRNA) = "RNA"
pdf("Brain_GSE188545_MTG_celltype_Marker1.pdf", width = 4.6,height = 6.7)
markers.to.plot <- c("PECAM1","VWF","FLT1","PLPP1","ARL15","MECOM","IGFBP3","VEGFC",
                     "RGCC","BTNL9","ITGA1", 
                     "ATP10A","SLC39A10","BBS9","IGF1R","ABCB1","FLI1","THSD4","GPCPD1",
                     "CLDN5","MT2A","BSG","ID1","CRIP2","LY6E","TAGLN2","BST2","S100A10") ##special capillary
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()

# PMID34582785_NAc
DefaultAssay(scRNA) = "RNA"
pdf("Brain_PMID34582785_NAc_celltype_Marker1.pdf", width = 4.2,height = 6.7)
markers.to.plot <- c("PECAM1","VWF","FLT1","PLPP1","ARL15","MECOM","IGFBP3","VEGFC",
                     "RGCC","BTNL9","ITGA1", 
                     "ATP10A","SLC39A10","BBS9","IGF1R","ABCB1","FLI1","THSD4","GPCPD1",
                     "CLDN5","MT2A","BSG","ID1","CRIP2","LY6E","TAGLN2","BST2","S100A10") ##special capillary
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()

# E_MTAB_8230
DefaultAssay(scRNA) = "RNA"
table(scRNA$Species)
scRNA1=scRNA[,scRNA$Species %in% c("macaque")]
pdf("Brain_E_MTAB_8230_PFC_macaque_celltype_Marker1.pdf", width = 4.2,height = 6.7)
markers.to.plot <- c("PECAM1","VWF","FLT1","PLPP1","ARL15","MECOM","IGFBP3","VEGFC",
                     "RGCC","BTNL9","ITGA1", 
                     "ATP10A","SLC39A10","BBS9","IGF1R","ABCB1","FLI1","THSD4","GPCPD1",
                     "CLDN5","MT2A","BSG","ID1","CRIP2","LY6E","TAGLN2","BST2","S100A10") ##special capillary
DotPlot(scRNA1, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()

# mouse
DefaultAssay(scRNA) = "SCT"
table(scRNA$Organ)
scRNA1=scRNA[,scRNA$Species %in% c("macaque")]
pdf("mouse_BSG_celltype_Marker1.pdf", width = 4.5,height = 5)
markers.to.plot <- c("Pecam1","Flt1","Plpp1","Arl15","Igfbp3",
                     "Rgcc","Btnl9","Itga1", 
                     "Atp10a","Slc39a10","Igf1r","Gpcpd1",
                     "Cldn5","Bsg","Id1","Crip2","Ly6e","Tagln2","Bst2","S100a10") ##special capillary
DotPlot(scRNA, features = markers.to.plot,group.by = "Organ") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()

# TS
scRNA=scRNA_predictions
DefaultAssay(scRNA) = "RNA"
table(scRNA$celltype)
pdf("TS_BSG_celltype_Marker1.pdf", width = 5.3,height = 3.5)
markers.to.plot <- c("CDH5","VWF","FLT1","ARL15","MECOM",
                     "RGCC", 
                     "ATP10A","IGF1R","GPCPD1",
                     "CLDN5","BSG","ID1","CRIP2","TAGLN2","BST2") ##special capillary
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()
################################################################################
# Extended data Fig. 7B
################################################################################
# GSM171668
pdf("Lung_GSM171668_celltype_Marker1.pdf", width = 4.5,height = 5.3)
markers.to.plot <- c("PECAM1","VWF","ARL15","VEGFC",
                     "ACKR1","CCL14","TLL1",
                     "RGCC","BTNL9","ITGA1","CCL21","MMRN1","LYVE1","TFPI","PROX1", #capillary
                     "HSPA1B","HSPA1A","F8","NR4A1", #artery
                     "AFF3","PRICKLE2","KHDRBS2","EDNRB","F2RL3") # capillary
DotPlot(scRNA, features = markers.to.plot,group.by = "manual_annotation_cluster") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()

################################################################################
# Extended data Fig. 7C
################################################################################
# GSE183852
load("/Users/niuruize/Downloads/scRNA/EC/classifier/heart/GSE183852/GSE183852_DCM_Integrated.Robj")
pdf("Heart_GSE183852_celltype_Marker1.pdf", width = 5,height = 7.5)
markers.to.plot <- c("PECAM1","VWF","PLPP1","ARL15","VEGFC","MECOM",
                     "ACKR1","CCL14",
                     "RGCC","BTNL9","ITGA1","CCL21","MMRN1","RELN","LYVE1","PROX1","TSPAN5",
                     "ID3","ADCY4", #artery
                     "PLXNA2","GFOD1","GRB10","F8","CYYR1","RAPGEF1","SEC14L1", #vein
                     "CAV1","IFITM3", #capillary
                     "SNTG2","CATSPERB","SMAD1","ADAM19","TSHZ2","TANC2","SEMA6A") #lymphatic
DotPlot(scRNA, features = markers.to.plot,group.by = "Names") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()
# zenodo6578553
scRNA <- readRDS("~/Downloads/scRNA/EC/classifier/heart/PMID35948637/MI.rds")
rownames(scRNA@assays$RNA@data) <- scRNA@assays$RNA@meta.features$feature_name
rownames(scRNA@assays$RNA@counts) <- scRNA@assays$RNA@meta.features$feature_name
DefaultAssay(scRNA) <- "RNA"
pdf("Heart_zenodo6578553_MI_celltype_Marker1.pdf", width = 4.2,height = 7.5)
markers.to.plot <- c("PECAM1","VWF","PLPP1","ARL15","MECOM","VEGFC","IGFBP3",
                     "ACKR1","CCL14","TLL1","FAM155A",
                     "RGCC","BTNL9","ITGA1",
                     "CCL21","MMRN1",
                     "LDB2","GMDS","ADCY4", #artery
                     "PLXNA2","GFOD1","GRB10","PLEKHG1","EVA1C","RAPGEF1","GALNT18","SEC14L1", #vein
                     "SNTG2","SMAD1","SEMA6A") #lymphatic
DotPlot(scRNA, features = markers.to.plot,group.by = "cell_type") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()

# zenodo6578553
scRNA <- readRDS("~/Downloads/scRNA/EC/classifier/heart/PMID35948637/normal.rds")
rownames(scRNA@assays$RNA@data)=scRNA@assays$RNA@meta.features$feature_name
rownames(scRNA@assays$RNA@counts)=scRNA@assays$RNA@meta.features$feature_name
pdf("Heart_zenodo6578553_Normal_celltype_Marker1.pdf", width = 4.6,height = 7.5)
markers.to.plot <- c("PECAM1","VWF","PLPP1","ARL15","MECOM","VEGFC","IGFBP3",
                     "ACKR1","CCL14","TLL1","FAM155A",
                     "RGCC","BTNL9","ITGA1",
                     "CCL21","MMRN1",
                     "LDB2","GMDS","ADCY4", #artery
                     "PLXNA2","GFOD1","GRB10","PLEKHG1","EVA1C","RAPGEF1","GALNT18","SEC14L1", #vein
                     "SNTG2","SMAD1","SEMA6A") #lymphatic
DotPlot(scRNA, features = markers.to.plot,group.by = "cell_type") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()
################################################################################
# Extended data Fig. 7D
################################################################################
# GSM4850577
scRNA <- scRNA[,scRNA$SampleID %in% c("Normal4")]
DefaultAssay(scRNA) <- "RNA"
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
pdf("Bladder_GSM4850577_celltype_Marker1.pdf", width = 4.1,height = 5.5)
markers.to.plot <- c("PECAM1","VWF","ARL15","VEGFC","MECOM",
                     "ACKR1","CCL14","ZNF385D","IL1R1",
                     "MMRN1","TFPI","PROX1",
                     "FLNB","CTTNBP2NL",
                     "CCL23","IFI27","IFITM3","SOCS3", #vein
                     "NR5A2","RASGRF2","PDE10A","CYYR1","SHANK3", #capillary
                     "HECW2") # lymphatic
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()
# GSE129845
scRNA <- scRNA[,scRNA$SampleID %in% c("Normal2","Normal3")]
DefaultAssay(scRNA) <- "RNA"
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
DefaultAssay(scRNA) = "RNA"
pdf("Bladder_GSE129845_celltype_Marker1.pdf", width = 4.1,height = 5.5)
markers.to.plot <- c("PECAM1","VWF","ARL15","VEGFC",
                     "ACKR1","CCL14","TLL1","CLU","ZNF385D",
                     "MMRN1","TFPI",
                     "EXOC6","CRIM1",
                     "CCL23","IFI27","CSF3","RNASE1","C2CD4B","SOCS3","TM4SF1", #vein
                     "TRIM5","NR5A2","RASGRF2","CYYR1", #capillary
                     "PLCG2")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()

# GSE135337
scRNA <- scRNA[,scRNA$SampleID %in% c("BCH1","BCH2","BCH3","BCH4","BCL1","BCL2","BCL3","BCN")]
DefaultAssay(scRNA) <- "RNA"
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
DefaultAssay(scRNA) = "RNA"
pdf("Bladder_GSE135337_celltype_Marker1.pdf",width = 4.1,height = 4.5)
markers.to.plot <- c("PECAM1","VWF","PLPP1","IGFBP3","MECOM",
                     "EXOC6","CRIM1",
                     "IFI27","IFITM3","B2M","BTF3","ITM2B","TM4SF1", #vein
                     "TRIM5","NR5A2","RASGRF2","CYYR1", #capillary
                     "PLCG2")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()

################################################################################
# Extended data Fig. 7E
################################################################################
# GSM171668
pdf("Liver_GSM171668_celltype_Marker1.pdf", width = 6.2,height = 5.7)
markers.to.plot <- c("PECAM1","VWF","PLPP1","ARL15","MECOM", #artery
                     "FAM155A","IL1R1", # vein
                     "BTNL9", #capillary
                     "CCL21","MMRN1","RELN","TFPI","COLEC12", #lymphatic
                     "BMPER","LDB2","WWTR1", #artery
                     "SMAD6","F8","GRB10","PLEKHG1", #vein
                     "SMAD1","NPAS3","SEMA6A") # lymphatic
DotPlot(scRNA, features = markers.to.plot,group.by = "predicted_celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()
pdf("Liver_GSM171668_celltype_Marker2.pdf", width = 11.2,height = 6.7)
markers.to.plot <- c("PECAM1","VWF","PLPP1","ARL15","MECOM", #artery
                     "FAM155A","IL1R1", # vein
                     "BTNL9", #capillary
                     "CCL21","MMRN1","RELN","TFPI","COLEC12", #lymphatic
                     "BMPER","LDB2","WWTR1", #artery
                     "SMAD6","F8","GRB10","PLEKHG1", #vein
                     "SMAD1","NPAS3","SEMA6A") # lymphatic
DotPlot(scRNA, features = markers.to.plot,group.by = "manual_annotation") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()
pdf("Liver_GSM171668_celltype_Marker3.pdf", width = 4.5,height = 4.7)
markers.to.plot <- c("PECAM1","VWF","PLPP1","ARL15","MECOM", #artery
                     "FAM155A", # vein
                     "BTNL9", #capillary
                     "CCL21","MMRN1","RELN","TFPI", #lymphatic
                     "BMPER","LDB2","WWTR1", #artery
                     "SMAD6","F8","GRB10","PLEKHG1", #vein
                     "NPAS3","SEMA6A") # lymphatic
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()

################################################################################
# Extended data Fig. 7F
################################################################################
# GSE132023
pdf("Kidney_GSE132023_celltype_Marker1.pdf", width = 3.5,height = 4.5)
markers.to.plot <- c("PECAM1","VWF","MECOM","ARL15", #artery
                     "RGCC","BTNL9","ITGA1", #capillary
                     "MMRN1","TFPI", #lymphatic
                     "FGF2", #artery
                     "KDR","FLT4","IGFBP5","EGLN1","ADM","CHST15", #vein
                     "TIMP3","CHRM3","RAPGEF3") # lymphatic
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()

# GSE109205
pdf("Kidney_GSE109205_celltype_Marker1.pdf", width = 4.3,height = 3.5)
markers.to.plot <- c("PECAM1","VWF","PLPP1", #artery
                     "RGCC","BTNL9", #capillary
                     "TFPI", #lymphatic"
                     "IGFBP5","KDR","FLT4","CEACAM1", #vein
                     "TIMP3","CHRM3") # lymphatic
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331")) + coord_flip()
dev.off()
