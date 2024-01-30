################################################################################
# Extended data Fig. 5
################################################################################

###--------------------------------------------Adipose--------------------------------------------###
scRNA=Adipose_Endo
##first
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=30, verbose=FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=20)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:20) %>% RunUMAP(reduction="harmony", dims = 1:20)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:20) %>% FindClusters(dims = 1:20, resolution = 0.1)
save(scRNA,file="Adipose_Harmony.RData")
##second
scRNA = scRNA[,scRNA$seurat_clusters %in% c("0","1","2","3","4","5")]
cellinfo <- subset(scRNA@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.rb","percent.HB",
                                                       "SampleID","Diagnosis","Age","Sex","Tangle.Stage","APOE","Organ","Region","Race","DoubletFinder","Group","datasets"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
##
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=20, verbose=FALSE)
ElbowPlot(scRNA, ndims = 20)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=20)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:20) %>% RunUMAP(reduction="harmony", dims = 1:20)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:20) %>% FindClusters(dims = 1:20, resolution = 0.2)
save(scRNA,file="Adipose_Harmony.RData")
##
current.cluster.ids <- c("0","1","2","3","4")
new.cluster.ids <- c("Ad2","Ad4","Ad1","Ad5","Ad3") 
scRNA$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype <- factor(scRNA$celltype,levels=c("Ad1","Ad2","Ad3","Ad4","Ad5"))
table(scRNA$seurat_clusters,scRNA$celltype)
new.cluster.ids <- c("capillary","vein","artery","lymphatic","cap.venous") 
scRNA$celltype_label <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label <- factor(scRNA$celltype_label,levels=c("artery","capillary","cap.venous","vein","lymphatic"))
table(scRNA$celltype_label,scRNA$celltype)
scRNA$celltype_label2 <- scRNA$celltype_label
save(scRNA,file="Adipose_Harmony.RData")
##
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
scRNA=Aorta_Endo
##first
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=30, verbose=FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=10)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:10) %>% RunUMAP(reduction="harmony", dims = 1:10)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:10) %>% FindClusters(dims = 1:10, resolution = 0.1)
save(scRNA,file="Aorta_Harmony.RData")
##second
table(scRNA$seurat_clusters)
scRNA = scRNA[,scRNA$seurat_clusters %in% c("1","2","3","4")]
cellinfo <- subset(scRNA@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.rb","percent.HB",
                                                       "SampleID","Diagnosis","Age","Sex","Tangle.Stage","APOE","Organ","Region","Race","DoubletFinder","Group","datasets"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=10, verbose=FALSE)
ElbowPlot(scRNA, ndims = 10)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=10)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:10) %>% RunUMAP(reduction="harmony", dims = 1:10)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:10) %>% FindClusters(dims = 1:10, resolution = 0.4)
save(scRNA,file="Aorta_Harmony.RData")
##
current.cluster.ids <- c("0","1","2")
new.cluster.ids <- c("Ao2","Ao1","Ao3") 
scRNA$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype <- factor(scRNA$celltype,levels=c("Ao1","Ao2","Ao3"))
table(scRNA$seurat_clusters,scRNA$celltype)
new.cluster.ids <- c("vein","artery","interferon") 
scRNA$celltype_label <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label <- factor(scRNA$celltype_label,levels=c("artery","vein","interferon"))
table(scRNA$celltype_label,scRNA$celltype)
scRNA$celltype_label2 <- scRNA$celltype_label
table(scRNA$celltype_label2,scRNA$celltype)
save(scRNA,file="Aorta_Harmony.RData")
##
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
scRNA=Bladder_Endo
## first
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=30, verbose=FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=10)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:10) %>% RunUMAP(reduction="harmony", dims = 1:10)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:10) %>% FindClusters(dims = 1:10, resolution = 0.2)
save(scRNA,file="Bladder_Harmony.RData")
##second
table(scRNA$seurat_clusters)
scRNA = scRNA[,scRNA$seurat_clusters %in% c("0","1","2","3","4","5")]
cellinfo <- subset(scRNA@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.rb","percent.HB",
                                                       "SampleID","Diagnosis","Age","Sex","Tangle.Stage","APOE","Organ","Region","Race","DoubletFinder","Group","datasets"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=10, verbose=FALSE)
ElbowPlot(scRNA, ndims = 10)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=10)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:10) %>% RunUMAP(reduction="harmony", dims = 1:10)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:10) %>% FindClusters(dims = 1:10, resolution = 0.1)
save(scRNA,file="Bladder_Harmony.RData")
##
current.cluster.ids <- c("31","35","47","48","NA")
new.cluster.ids <- c("Adult","Adult","Aging","Aging","NA") 
scRNA$Group <- plyr::mapvalues(x = as.character(scRNA@meta.data$Age), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA$Group)
current.cluster.ids <- c("0","1","2","3","4","5")
new.cluster.ids <- c("Bl5","Bl2","Bl4","Bl6","Bl3","Bl1") 
scRNA$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype <- factor(scRNA$celltype,levels=c("Bl1","Bl2","Bl3","Bl4","Bl5","Bl6"))
table(scRNA$seurat_clusters,scRNA$celltype)
new.cluster.ids <- c("vein2","cap.arterial","vein1","lymphatic","capillary","artery") 
scRNA$celltype_label <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label <- factor(scRNA$celltype_label,levels=c("artery","cap.arterial","capillary","vein1","vein2","lymphatic"))
table(scRNA$celltype_label,scRNA$celltype)
new.cluster.ids <- c("vein","cap.arterial","vein","lymphatic","capillary","artery") 
scRNA$celltype_label2 <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label2 <- factor(scRNA$celltype_label2,levels=c("artery","cap.arterial","capillary","vein","lymphatic"))
table(scRNA$celltype_label2,scRNA$celltype)
save(scRNA,file="Bladder_Harmony.RData")
##
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
scRNA=Brain_Endo
scRNA = scRNA[,scRNA$Region %in% c("ACC","EC","HIP","MB","OC","OTC","PFC")]
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=20, verbose=FALSE)
ElbowPlot(scRNA, ndims = 20)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=10)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:10) %>% RunUMAP(reduction="harmony", dims = 1:10)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:10) %>% FindClusters(dims = 1:10, resolution = 0.1)
save(scRNA,file="Brain_Harmony.RData")
##second
table(scRNA$seurat_clusters)
scRNA = scRNA[,scRNA$seurat_clusters %in% c("0","1","2","3","4","5")]
cellinfo <- subset(scRNA@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.rb","percent.HB",
                                                       "SampleID","Diagnosis","Age","Sex","Tangle.Stage","APOE","Organ","Region","Race","Group","datasets"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=10, verbose=FALSE)
ElbowPlot(scRNA, ndims = 10)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=10)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:10) %>% RunUMAP(reduction="harmony", dims = 1:10)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:10) %>% FindClusters(dims = 1:10, resolution = 0.1)
save(scRNA,file="Brain_Harmony2.RData")
##
current.cluster.ids <- c("0","1","2","3","4","5")
new.cluster.ids <- c("Br5","Br2","Br3","Br4","Br1","Br6") 
scRNA$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype <- factor(scRNA$celltype,levels=c("Br1","Br2","Br3","Br4","Br5","Br6"))
table(scRNA$seurat_clusters,scRNA$celltype)
new.cluster.ids <- c("vein3","capillary","vein1","vein2","artery","RGS5+vein") 
scRNA$celltype_label <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label <- factor(scRNA$celltype_label,levels=c("artery","capillary","vein1","vein2","vein3","RGS5+vein"))
table(scRNA$celltype_label,scRNA$celltype)
new.cluster.ids <- c("vein","capillary","vein","vein","artery","RGS5+vein") 
scRNA$celltype_label2 <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label2 <- factor(scRNA$celltype_label2,levels=c("artery","capillary","vein","RGS5+vein"))
table(scRNA$celltype_label2,scRNA$celltype)
save(scRNA,file="Brain_Harmony2.RData")
##
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
scRNA=SC_Endo
scRNA = scRNA[,scRNA$Region %in% c("SC")]
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=20, verbose=FALSE)
ElbowPlot(scRNA, ndims = 20)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=10)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:10) %>% RunUMAP(reduction="harmony", dims = 1:10)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:10) %>% FindClusters(dims = 1:10, resolution = 0.1)
save(scRNA,file="SC_Harmony.RData")
##
current.cluster.ids <- c("0","1","2","3","4")
new.cluster.ids <- c("Sc3","Sc1","Sc4","Sc2","Sc5") 
scRNA$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype <- factor(scRNA$celltype,levels=c("Sc1","Sc2","Sc3","Sc4","Sc5"))
table(scRNA$seurat_clusters,scRNA$celltype)
new.cluster.ids <- c("capillary2","artery","vein1","capillary1","vein2") 
scRNA$celltype_label <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label <- factor(scRNA$celltype_label,levels=c("artery","capillary1","capillary2","vein1","vein2"))
table(scRNA$celltype_label,scRNA$celltype)
new.cluster.ids <- c("capillary","artery","vein","capillary","vein") 
scRNA$celltype_label2 <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label2 <- factor(scRNA$celltype_label2,levels=c("artery","capillary","vein"))
table(scRNA$celltype_label2,scRNA$celltype)
save(scRNA,file="SC_Harmony.RData")
##
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
scRNA=Esophagus_Endo
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims = 1:10)
scRNA <- FindClusters(scRNA, resolution = 0.2)
scRNA <- RunUMAP(scRNA, dims = 1:10)
save(scRNA,file="Esophagus_Harmony.RData")
##
current.cluster.ids <- c("0","1","2","3","4","5")
new.cluster.ids <- c("E6","E4","E5","E3","E1","E2") 
scRNA$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype <- factor(scRNA$celltype,levels=c("E1","E2","E3","E4","E5","E6"))
table(scRNA$seurat_clusters,scRNA$celltype)
new.cluster.ids <- c("vein2","capillary2","vein1","capillary1","artery1","artery2") 
scRNA$celltype_label <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label <- factor(scRNA$celltype_label,levels=c("artery1","artery2","capillary1","capillary2","vein1","vein2"))
table(scRNA$celltype_label,scRNA$celltype)
new.cluster.ids <- c("vein","capillary","vein","capillary","artery","artery") 
scRNA$celltype_label2 <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label2 <- factor(scRNA$celltype_label2,levels=c("artery","capillary","vein"))
table(scRNA$celltype_label2,scRNA$celltype)
save(scRNA,file="Esophagus_Harmony.RData")
##
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
scRNA=Heart_Endo
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=30, verbose=FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=20)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:20) %>% RunUMAP(reduction="harmony", dims = 1:20)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:20) %>% FindClusters(dims = 1:20, resolution = 0.1)
save(scRNA,file="Heart_Harmony.RData")
##second
scRNA = scRNA[,scRNA$seurat_clusters %in% c("0","1","2","3","4","5","6","8")]
cellinfo <- subset(scRNA@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.rb","percent.HB",
                                                       "SampleID","Diagnosis","Age","Sex","Tangle.Stage","APOE","Organ","Region","Race","DoubletFinder","Group","datasets"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=20, verbose=FALSE)
ElbowPlot(scRNA, ndims = 20)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=10)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:10) %>% RunUMAP(reduction="harmony", dims = 1:10)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:10) %>% FindClusters(dims = 1:10, resolution = 0.1)
save(scRNA,file="Heart_Harmony.RData")
current.cluster.ids <- c("0","1","2","3","4","5")
new.cluster.ids <- c("He2","He5","He3","He1","He4","He6") 
scRNA$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype <- factor(scRNA$celltype,levels=c("He1","He2","He3","He4","He5","He6"))
table(scRNA$seurat_clusters,scRNA$celltype)
new.cluster.ids <- c("cap.arterial1","vein","cap.arterial2","artery","capillary","lymphatic") 
scRNA$celltype_label <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label <- factor(scRNA$celltype_label,levels=c("artery","cap.arterial1","cap.arterial2","capillary","vein","lymphatic"))
table(scRNA$celltype_label,scRNA$celltype)
new.cluster.ids <- c("cap.arterial","vein","cap.arterial","artery","capillary","lymphatic") 
scRNA$celltype_label2 <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label2 <- factor(scRNA$celltype_label2,levels=c("artery","cap.arterial","capillary","vein","lymphatic"))
table(scRNA$celltype_label2,scRNA$celltype)
save(scRNA,file="Heart_Harmony.RData")
##
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
scRNA=MG_Endo
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=30, verbose=FALSE)
ElbowPlot(scRNA, ndims = 30)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=10)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:10) %>% RunUMAP(reduction="harmony", dims = 1:10)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:10) %>% FindClusters(dims = 1:10, resolution = 0.2)
save(scRNA,file="MG_Harmony.RData")
##second
scRNA = scRNA[,scRNA$seurat_clusters %in% c("0","2","3","4","6")]
cellinfo <- subset(scRNA@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.rb","percent.HB",
                                                       "SampleID","Diagnosis","Age","Sex","Tangle.Stage","APOE","Organ","Region","Race","DoubletFinder","Group","datasets"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
scRNA <- SCTransform(scRNA)
scRNA <- RunPCA(scRNA, npcs=20, verbose=FALSE)
ElbowPlot(scRNA, ndims = 20)
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=10)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:10) %>% RunUMAP(reduction="harmony", dims = 1:10)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:10) %>% FindClusters(dims = 1:10, resolution = 0.2)
save(scRNA,file="MG_Harmony.RData")
current.cluster.ids <- c("0","1","2","4")
new.cluster.ids <- c("MG3","MG4","MG2","MG1") 
scRNA$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype <- factor(scRNA$celltype,levels=c("MG1","MG2","MG3","MG4"))
table(scRNA$seurat_clusters,scRNA$celltype)
new.cluster.ids <- c("vein","interferon","capillary","artery") 
scRNA$celltype_label <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA$celltype_label <- factor(scRNA$celltype_label,levels=c("artery","capillary","vein","interferon"))
table(scRNA$celltype_label,scRNA$celltype)
scRNA$celltype_label2 <- scRNA$celltype_label
save(scRNA,file="MG_Harmony.RData")
##

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

##
pdf("Uterus_celltype_Marker1.pdf", width = 5.9,height = 2.4)
markers.to.plot <- c("FBLN5","GJA5","BMX","EFNB2","ADGRF5","EMCN","SGK1","RGCC","CA4","PRX","ACKR1","NR2F2","EPHB4","VCAM1","IFITM1","IFITM2","IFITM3","MKI67","TOP2A","CENPF")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()







