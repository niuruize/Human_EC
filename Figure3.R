################################################################################
# Figure 3
################################################################################
# Load R libraries
library(Seurat)
 ###--------------------------------------------Adipose--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_group.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Adipose_celltype_Marker1.pdf", width = 4.7,height = 1.6)
markers.to.plot <- c("FBLN5","GJA5","BMX","EFNB2","RGCC","CA4","SPARC","NR2F2","EPHB4","ACKR1","PROX1","LYVE1","IFITM1","IFITM2","IFITM3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###-------------------------------------------Aorta--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Aorta_celltype_Marker1.pdf", width = 4.5,height = 1.3)
markers.to.plot <- c("FBLN5","GJA5","BMX","EFNB2","VWF","EPHB4","EMCN","IFITM1","IFITM2","IFITM3","IFIT1","IFIT2","IFIT3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###--------------------------------------------Bladder--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Bladder_celltype_Marker1.pdf", width = 5.4,height = 2)
markers.to.plot <- c("FBLN5","GJA5","BMX","EFNB2","RGCC","CA4","EMCN","ADGRF5","VWF","NR2F2","EPHB4","ICAM1","ACKR1","PROX1","LYVE1","IFITM1","IFITM2","IFITM3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###-------------------------------------------Brain--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Brain_celltype_Marker1.pdf", width = 4.5,height = 2)
markers.to.plot <- c("FBLN5","GJA5","BMX","EFNB2","RGCC","CA4","VWF","EPHB4","ZNF385D","EMCN","EBF1","RGS5","IFITM1","IFITM2","IFITM3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###-------------------------------------------SC--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=1, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=1, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("SC_celltype_Marker1.pdf", width = 4.7,height = 2)
markers.to.plot <- c("FBLN5","GJA5","BMX","EFNB2","RGCC","CA4","EPHB4","EMCN","VWF","ZNF385D","IFITM1","IFITM2","IFITM3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###--------------------------------------------Esophagus--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size= 1, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=1, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Eso_celltype_Marker1.pdf", width = 5.2,height = 2.1)
markers.to.plot <- c("FBLN5","GJA5","EFNB2","BMX","LTBP4","LTBP1","NOTCH2","RGCC","CA4","NR2F2","EPHB4","ICAM1","ACKR1","IFITM1","IFITM2","IFITM3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###--------------------------------------------Heart--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Heart_celltype_Marker1.pdf", width = 5.2,height = 2.1)
markers.to.plot <- c("BMX","HES4","ADGRF5","ITGA1","DLL4","EFNB2","RGCC","CA4","VWF","EPHB4","PROX1","LYVE1","IFITM1","IFITM2","IFITM3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###------------------------------------------------MG--------------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("MG_celltype_Marker1.pdf", width = 5,height = 1.5)
markers.to.plot <- c("FBLN5","GJA5","EFNB2","BMX","RGCC","SPARC","CA4","NR2F2","EPHB4","ICAM1","ACKR1","IFITM1","IFITM2","IFITM3","IFIT1","IFIT2","IFIT3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###--------------------------------------------------Muscle--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Muscle_celltype_Marker1.pdf", width = 4.9,height = 2)
markers.to.plot <- c("FBLN5","GJA4","EFNB2","BMX","EMCN","RGCC","CA4","NR2F2","EPHB4","ICAM1","ACKR1","PROX1","LYVE1")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###--------------------------------------------------------GI--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size= 2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size= 2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("GI_celltype_Marker1.pdf", width = 4.2,height = 1.3)
markers.to.plot <- c("FBLN5","GJA5","EFNB2","DLL4","VWF","EPHB4","LYVE1","PROX1","IFITM1","IFITM2","IFITM3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###------------------------------------------------------Kidney--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Kidney_celltype_Marker1.pdf", width = 4.1,height = 1.5)
markers.to.plot <- c("EFNB2","DLL4","FBLN5","GJA5","SGK1","AFF3","NR2F2","IGFBP7","PLAT","TSPAN7","IFITM1","IFITM2","IFITM3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###-------------------------------------------------------Liver--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Liver_celltype_Marker1.pdf", width = 4.9,height = 1.8)
markers.to.plot <- c("EFNB2","DLL4","GJA4","CLEC4G","CLEC4M","FCN2","VWF","NR2F2","EPHB4","IFITM1","IFITM2","IFITM3","RELN","STAB2")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###-------------------------------------------Lung--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Lung_celltype_Marker1.pdf", width = 4.5,height = 2)
markers.to.plot <- c("EFNB2","DLL4","BMX","SPARC","RGCC","EMCN","CA4","VWF","NR2F2","EPHB4","IFITM1","IFITM2","IFITM3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###-------------------------------------------Decidua--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=1, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=1, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Decidua_celltype_Marker1.pdf", width = 4.4,height = 1.5)
markers.to.plot <- c("EFNB2","DLL4","BMX","RGCC","SPARC","EPHB4","VCAM1","IFITM1","IFITM2","IFITM3","PROX1","LYVE1")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###-------------------------------------------Placenta--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Placenta_celltype_Marker1.pdf", width = 4.4,height = 1.8)
markers.to.plot <- c("GJA4","HEY1","NR2F2","EPHB4","ICAM1","ACKR1","SELP","IFITM1","IFITM2","IFITM3","LYVE1")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###-------------------------------------------Ovary--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Ovary_celltype_Marker1.pdf", width = 5.3,height = 2.1)
markers.to.plot <- c("EFNB2","DLL4","BMX","CA4","SPARC","ADGRL2","ICAM1","VCAM1","ACKR1","IFITM1","IFITM2","IFITM3","IFIT1","IFIT2","IFIT3","PROX1","LYVE1")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###--------------------------------------------Pancrease--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=1, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=1, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Pancreas_celltype_Marker1.pdf", width = 4.7,height = 1.7)
markers.to.plot <- c("EFNB2","FBLN5","BMX","RGCC","CA4","SPARC","SGK1","EMCN","EPHB4","NR2F2","VCAM1","ICAM1","PROX1","LYVE1")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###--------------------------------------------PIT--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("PIT_celltype_Marker1.pdf", width = 4.6,height = 1.5)
markers.to.plot <- c("EFNB2","CLU","MECOM","BMX","NOTCH1","RGCC","CA4","EPHB4","IFITM1","IFITM2","IFITM3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###--------------------------------------------Retina--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Retina_celltype_Marker1.pdf", width = 4.7,height = 1.2)
markers.to.plot <- c("NOTCH1","EFNB2","FBLN5","DLL4","RGCC","NR2F2","EPHB4","IFITM1","IFITM2","IFITM3","IFIT1","IFIT2","IFIT3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###--------------------------------------------Skin--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Skin_celltype_Marker1.pdf", width = 5.5,height = 2.3)
markers.to.plot <- c("NOTCH1","ELN","MECOM","EFNB2","DLL4","FBLN5","IGFBP3","ARL15","RGCC","NR2F2","EPHB4","IFITM1","IFITM2","IFITM3","PROX1","LYVE1")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###--------------------------------------------Testis--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Testis_celltype_Marker1.pdf", width = 4.8,height = 1.7)
markers.to.plot <- c("FBLN5","GJA5","GJA4","VWF","EMCN","SELP","IFITM1","IFITM2","IFITM3","IFIT1","IFIT2","IFIT3","PROX1","LYVE1")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()


###--------------------------------------------Thyroid--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Thyroid_celltype_Marker1.pdf", width = 5.0,height = 1.7)
markers.to.plot <- c("FBLN5","GJA5","GJA4","EFNB2","EMCN","RGCC","CA4","VWF","NR2F2","EPHB4","PROX1","LYVE1","IFITM1","IFITM2","IFITM3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###--------------------------------------------Trachea--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Trachea_celltype_Marker1.pdf", width = 5.0,height = 1.6)
markers.to.plot <- c("SOX17","GJA4","DLL4","EFNB2","EMCN","RGCC","EMP2","CA4","EPHB4","NR2F2","ICAM1","IFITM1","IFITM2","IFITM3")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

###--------------------------------------------Uterus--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_label", pt.size=0.2, label = TRUE,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_label.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')

pdf("Uterus_celltype_Marker1.pdf", width = 5.9,height = 2.4)
markers.to.plot <- c("FBLN5","GJA5","BMX","EFNB2","ADGRF5","EMCN","SGK1","RGCC","CA4","PRX","ACKR1","NR2F2","EPHB4","VCAM1","IFITM1","IFITM2","IFITM3","MKI67","TOP2A","CENPF")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()
