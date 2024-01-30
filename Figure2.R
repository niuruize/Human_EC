################################################################################
#  pySCENIC
################################################################################
#=====================================================================
#                        1.Matrix of expression
#=====================================================================
library(Seurat)
dim(scRNA)
table(Idents(scRNA))
Idents(scRNA) <- "Organ"
table(scRNA$Organ)
#scRNA = subset(scRNA, idents = "Low_quality", invert = TRUE)
scRNA$Organ <- as.factor(as.character(scRNA$Organ))
table(scRNA$Organ)
#scRNA$Group_Organ <- paste(scRNA$Group, scRNA$Organ, sep = "_")
#Idents(scRNA) <- 'Group_Organ'
table(Idents(scRNA))
scRNA <- subset(scRNA, subset = nFeature_RNA > 500)
table(Idents(scRNA))
sce_test <- subset(scRNA, downsample = 500)
table(Idents(sce_test))
table(sce_test$Organ)
dim(sce_test)
Idents(sce_test) <- 'Organ'
write.csv(t(as.matrix(sce_test@assays$SCT@data)),file = "sce_exp.csv")
colnames(sce_test@meta.data)
cellInfo <- sce_test@meta.data
colnames(cellInfo) <- c('Organ','Region','Age','Group','Sex','copykat.pred')
head(cellInfo)
write.csv(cellInfo, file = "cellInfo.csv")
saveRDS(sce_test, file = 'sce_test.rds')
###
library(Seurat)
library(Matrix)
seurat_to_adata <- function(object,#seurat对象
                            Dimension=c('UMAP','TSNE'),#降维方式
                            path){#文件保存路径
  seurat_obj <- object
  seurat_obj$barcode <- colnames(seurat_obj)
  if(Dimension=='UMAP'){
    cell.embeddings<- seurat_obj@reductions$umap@cell.embeddings
    seurat_obj$UMAP_1 <- cell.embeddings[,1]
    seurat_obj$UMAP_2 <- cell.embeddings[,2]
  }else{
    
    cell.embeddings<- seurat_obj@reductions$tsne@cell.embeddings
    seurat_obj$TSNE_1 <- cell.embeddings[,1]
    seurat_obj$TSNE_2 <- cell.embeddings[,2]
  }
  #metadat
  write.csv(seurat_obj@meta.data, file=paste0(path,'metadata.csv'), quote=F, row.names=F)
  #matrix
  counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
  writeMM(counts_matrix, file=paste0(path, 'counts.mtx'))
  #PCA
  write.csv(seurat_obj@reductions$pca@cell.embeddings, file=paste0(path,'pca.csv'), quote=F,row.names=F)
  
  #gene name
  write.table(data.frame('gene'=rownames(counts_matrix)),file=paste0(path,'gene_names.csv'),
              quote=F,row.names=F,col.names=F)
}
#
seurat_to_adata(sce_test,Dimension='TSNE',path = './seurat_')
#=====================================================================
#                        2.matrix to loom
#=====================================================================
vim trans.py
import os, sys
os.getcwd()
os.listdir(os.getcwd())
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("sce_exp.csv")
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("sce.loom",x.X.transpose(),row_attrs,col_attrs)
:wq
python trans.py
ls
#=====================================================================
#                        3.pyscenic
#=====================================================================
nohup pyscenic grn --num_workers 16 \
                --sparse \
                --method grnboost2 \
                --output grn.csv \
                sce.loom \
                allTFs_hg38.txt &
pyscenic grn --num_workers 8 \
--sparse \
--method grnboost2 \
--output grn.csv \
sce.loom \
/Users/niuruize/Downloads/scRNA/pyscenic/allTFs_hg38.txt
#
#pip3 uninstall numpy
#pip3 install numpy==1.21.1 -i https://pypi.tuna.tsinghua.edu.cn/simple --default-timeout=600 --force-reinstall 
#conda install numpy==1.21.1 --force-reinstall 
nohup pyscenic ctx --num_workers 16 \
  --output sce.regulons.csv \
  --expression_mtx_fname sce.loom \
  --all_modules \
  --mask_dropouts \
  --mode "dask_multiprocessing" \
  --min_genes 10 \
  --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  grn.csv \
  hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather &
#
nohup pyscenic aucell --num_workers 16 \
  --output sce_SCENIC.loom \
  sce.loom \
  sce.regulons.csv &
#=====================================================================
#                        4.visualization
#=====================================================================
# Figure 2A
# The code for the visualization is shown in pyscenic.ipynb
# Figure 2B module regulon
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
sce_SCENIC <- open_loom("sce_SCENIC.loom")
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
data <- regulonAUC@assays@data$AUC
data <- t(data)
#rm(sce_SCENIC, regulonAUC, data)
regulon <- data
# calculate CSI
CSI_matrix_cal <- function(regulon,
                           CSI_threshold,
                           module_k,
                           module_color=F,
                           Heatmap_col=NULL,
                           legend_parm = c("number","character"),
                           rect_color,
                           label_reg=NULL){
  
  #calculate CSI
  
  Mcor<-cor(regulon)
  n<-nrow(Mcor)

  CSI<-matrix(nrow=n,ncol=n)
  
  for (i in 1:n){
    for(j in 1:n){
      
      if(i==j) {
        
        CSI[i,j] <- 1
        
      } else{
        
        nodeA <- names(which(Mcor[i,]>= Mcor[i,j]-0.05))
        nodeB <- names(which(Mcor[,j]>= Mcor[i,j]-0.05))
        CSI[i,j]<- 1-((length(unique(c(nodeA,nodeB))))/n)
          
      }
      
    }
    
  }
  
  rownames(CSI)<-colnames(regulon)
  colnames(CSI)<-colnames(regulon)
  
  CSI_matrix <- as.data.frame(CSI)
  write.csv(CSI_matrix, file = './CSI_matrix.csv')
  
  
  #Heatmap-draw
  require(ComplexHeatmap)
  require(pheatmap)
  
  CSI[CSI <= CSI_threshold]=0
  
  
  if(is.null(Heatmap_col)){
    
    col<-colorRampPalette(c("#FAF9DA","#28245F"))(100)
    
  }else{
    
    col = Heatmap_col
    
  }
  
  
  x=pheatmap::pheatmap(CSI,
                       color=col,
                       clustering_method = "ward.D2",
                       show_rownames=FALSE,
                       show_colnames = FALSE,
                       cutree_rows = module_k,
                       cutree_cols = module_k)
  
  
  annotation_row <- data.frame(Cluster=factor(cutree(x$tree_row, module_k)))
  annotation_col <- data.frame(Cluster=factor(cutree(x$tree_col, module_k)))
  
  row_order <- annotation_row
  row_order$regulon <- rownames(row_order)
  row_order <- row_order[order(row_order$Cluster),]
  write.csv(row_order, file = "./Module.csv")
  
  
  anno_col = annotation_col
  anno_col$TF <- rownames(anno_col)
  
  index <- x$tree_col$order
  TFs <- x$tree_col$labels
  ord_TF <- c()
  for (i in index) {
    
    ord_TF <- append(ord_TF, TFs[i])
    
  }
  
  anno_col <- anno_col[ord_TF,]
  anno_col$Modules <- paste0("Module",anno_col$Cluster)
  
  
  
  if(module_color==F){
    
    calm = c("#7DD06F", "#844081", "#688EC1", "#C17E73", "#484125", 
            "#6CD3A7", "#597873","#7B6FD0", "#CF4A31", "#D0CD47",
            "#722A2D", "#CBC594", "#D19EC4", "#5A7E36", "#D4477D",
            "#403552", "#76D73C", "#96CED5", "#CE54D1", "#C48736")
                      
    module_num <- unique(anno_col$Cluster) 
    cluster_color = setNames(calm[1:module_k],module_num) 
    
    
    
  }else{
    
    
    module_num <- unique(anno_col$Cluster)
    cluster_color = setNames(module_color,module_num) 
    
    
  }
  
  
  cluster_color_m <- as.data.frame(cluster_color)
  cluster_color_m$Modules <- paste0("Module",rownames(cluster_color_m))
  rownames(cluster_color_m) <- cluster_color_m$Modules
  

  if(legend_parm == "number"){
    
    heatmap_legend_param = list(color_bar = "continuous",
                                legend_direction = "vertical",
                                legend_width = unit(1, "cm"),
                                legend_height = unit(5, "cm"),
                                title = "Connection specificity index (CSI)",
                                title_position="leftcenter-rot",
                                border ="black",
                                at = c(0,0.2,0.4,0.6,0.8,1),
                                labels = c(0,0.2,0.4,0.6,0.8,1),
                                labels_gp = gpar(fontsize = 8,col='black',font = 3))
    
  }
  
  
  if(legend_parm == "character"){
    
    heatmap_legend_param = list(color_bar = "continuous",
                                legend_direction = "vertical",
                                legend_width = unit(1, "cm"),
                                legend_height = unit(5, "cm"),
                                title = "Connection specificity index (CSI)",
                                title_position="leftcenter-rot",
                                border ="black",
                                at = c(0,0.5,1),
                                labels = c("low","mid","high"),
                                labels_gp = gpar(fontsize = 8,col='black',font = 3))
  }
  
  hm = ComplexHeatmap::pheatmap(CSI, #ph$tree_row$order
                          annotation_row=annotation_row,
                          annotation_col=annotation_col,
                          clustering_method = "ward.D2",
                          show_rownames=FALSE,
                          show_colnames = FALSE,
                          color=col,
                          name = "ht",
                          treeheight_row = 20,
                          treeheight_col = 20,
                          annotation_names_col = F,
                          annotation_names_row = F,
                          annotation_legend=F,
                          annotation_colors= list(Cluster = cluster_color),
                          heatmap_legend_param = heatmap_legend_param)
  
  
  if(!is.null(label_reg)){
    
   label_reg <- label_reg
   index <- which(rownames(CSI)%in% label_reg)
   customRowLabel <- rownames(CSI)[index]

   hm <- hm+rowAnnotation(anno = anno_mark(at = index,
                                 labels = customRowLabel,
                                 side = "right",
                                 padding = unit(2, "mm"),
                                 link_width = unit(3, "mm"),
                                 extend = unit(0.1, "mm"),
                                 labels_gp = gpar(fontsize = 8),
                                 link_gp = gpar(col='black',lwd=1)))
    
  }
  
  
  draw(hm)
  
  ord = anno_col$Cluster
  dup = (which(!duplicated(ord)) - 1)
  fract = dup / nrow(anno_col)
  width =  c(fract[-1], 1) - fract
  
  decorate_heatmap_body("ht", {
    grid.rect(unit(fract, "native"), 
              unit(1-fract, "native"), 
              unit(width, "native"), 
              unit(width, "native"), 
              hjust = 0, 
              vjust = 1, 
              gp = gpar(col = rect_color, lty = 1, lwd = 2, fill=NA))
  })
  
  
  label_m <- unique(anno_col$Modules)
  cluster_color_m <- cluster_color_m[label_m, ]
  
  decorate_heatmap_body("ht", {
    
    
    grid.text(label_m, 
              unit(fract+0.25, "native"), 
              unit(1-fract-0.05, "native"), 
              gp=gpar(fontsize=15, col=cluster_color_m$cluster_color, fontface="bold"))
    
  })
  
  
  return(hm)
}

# module_TF <- CSI_matrix_cal(regulon = regulon,
#                             CSI_threshold = 0.5,
#                             module_k = 6,
#                             legend_parm = "number",
#                             rect_color="red",
#                             label_reg = label_reg)
# module_TF

pdf('module3.pdf', width=8, height=6.5)
module_TF <- CSI_matrix_cal(regulon = regulon,
                            CSI_threshold = 0.5,
                            module_k = 6,
                            legend_parm = "number",
                            rect_color="red",
                            label_reg = label_reg$label_reg)
#draw(module_TF)
dev.off()

# grid.text("Module2", x=0.34, y=0.8, gp=gpar(fontsize=15, col="#844081", fontface="bold"))
# grid.text("Module1", x=0.5, y=0.6, gp=gpar(fontsize=15, col="#7DD06F", fontface="bold"))
# grid.text("Module3", x=0.65, y=0.45, gp=gpar(fontsize=15, col="#688EC1", fontface="bold"))
# grid.text("Module4", x=0.75, y=0.3, gp=gpar(fontsize=15, col="#C17E73", fontface="bold"))
# grid.text("Module6", x=0.85, y=0.2, gp=gpar(fontsize=15, col="#6CD3A7", fontface="bold"))
# 


# Mcor<-cor(regulon)
# Mcor[Mcor <= 0.5]=0
# col<-colorRampPalette(c("#FAF9DA","#28245F"))(100)
# pheatmap::pheatmap(Mcor,
#                    color=col,
#                    clustering_method = "ward.D2",
#                    show_rownames=FALSE,
#                    show_colnames = FALSE,
#                    cutree_rows = 6,
#                    cutree_cols = 6)
#label regulon in different module
module <- read.csv("Module.csv", header = T, row.names = 1,check.names = F)
label_reg <- c("TCF3(-)","PRDM6(-)","NFKB1(+)","SOX4(+)","PGR(-)","ZIC3(+)",
               "HOXD10(-)","HOXD9(+)","BCL3(+)","KLF4(+)","ZNF777(-)","EHF(-)",
                            "IRF8(+)","PBX1(+)","IRX5(-)","HNF4A(+)","GATA4(+)",
               "FOXF1(-)","FOXF1(+)","HMGA1(+)","ILF2(-)","SOX7(+)","PPARG(+)",
               "ZNF135(-)","MAFB(-)","HNF4G(-)",                    "ZBTB18(+)",
               "ARID3A(+)","HOXA10(-)","DDIT3(-)", "ZNF253(-)",      
               "MYRF(-)",             "HOXA7(-)",            "WT1(-)","HOXD9(-)",
               "HES5(+)","NKX2-1(-)","LTF(+)",                "FOS(+)","OSR2(+)") 

label_reg <- as.data.frame(label_reg)
CSI <- read.csv("CSI_matrix.csv", header = T,row.names = 1,check.names = F)
pdf('module2.pdf', width=8, height=6)
hm_tf <- module_TF+rowAnnotation(link = anno_mark(at = which(rownames(CSI) %in% label_reg$label_reg), 
                                         labels = label_reg$label_reg, labels_gp = gpar(fontsize = 8)))

hm_tf
draw(hm_tf) 
dev.off()
#
CSI <- read.csv("CSI_matrix.csv", header = T,row.names = 1,check.names = F)
col<-colorRampPalette(c("#FAF9DA","#28245F"))(100)
CSI[CSI <= 0.5]=0
CSI <- as.matrix(CSI)


ComplexHeatmap::pheatmap(CSI, 
                         clustering_method = "ward.D2",
                         show_rownames=FALSE,
                         show_colnames = FALSE,
                         color=col,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         annotation_legend=F,
                         cutree_rows = 6,
                         cutree_cols = 6,
                         layer_fun = function(j, i, x, y, width, height, fill) {
                           v = pindex(CSI, i, j)
                           if(sum(v > 0)/length(v) > 0.4) {
                             grid.rect(gp = gpar(lwd = 3, fill = "transparent", col='orange'))
                           }
                         })
# Figure 2B wordcloud
# Organ-regulon-rank
library(wordcloud2);library(webshot);library(htmlwidgets);library(Ipaper)
module <- read.csv("Module.csv", header = T, row.names = 1, check.names = F)

#module regulon
module_reg <- list()
for (i in 1:6) {
  module_r <- subset(module, Cluster==i)
  module_r <- module_r$regulon
  
  module_reg[[i]] <- module_r
}
#Module-auc
Average_mr <- matrix(nrow=nrow(regulon),ncol=6)
rownames(Average_mr) <- rownames(regulon)

for (i in 1:6) {
  
  Average_M <- regulon[,module_reg[[i]]]
  Average_M$module_A <- rowMeans(Average_M)
  Average_mr[,i] <- Average_M$module_A
  
}
colnames(Average_mr) <- c("M1","M2","M3","M4","M5","M6")
library(Seurat)
sce_test <- readRDS("sce_test.rds")
reduc.df <- as.data.frame(Embeddings(sce_test$Organ))
reduc.df <- cbind(reduc.df, Average_mr)

Module_mean <- aggregate(reduc.df[,2:7], by=list(type=reduc.df$`sce_test@meta.data$Organ`),mean)
wordcloud2(Module_mean[,c(1,7)],shape = 'cardioid',size=0.5) + WCtheme(2)

color = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
          '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
          '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
          '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
          '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
          '#968175','#985a38', '#aac6d2')

#install.packages("ggwordcloud")
library(ggwordcloud)

#geom_text_wordcloud_ares
wordf <- arrange(Module_mean[,c(1,2)], -M1)
p1=ggplot(wordf,aes(label=type,size=M1,color=color))+
  geom_text_wordcloud_area(shape="circle")+
  scale_size_area(max_size=20)+
  theme_minimal()
ggsave(filename = "regulon_M1.pdf", plot = p1, device = 'pdf', width = 14, height = 10, units = 'cm')
#
wordf <- arrange(Module_mean[,c(1,3)], -M2)
p1=ggplot(wordf,aes(label=type,size=M2,color=color))+
  geom_text_wordcloud_area(shape="circle")+
  scale_size_area(max_size=20)+
  theme_minimal()
ggsave(filename = "regulon_M2.pdf", plot = p1, device = 'pdf', width = 14, height = 10, units = 'cm')
#
wordf <- arrange(Module_mean[,c(1,4)], -M3) 
p1=ggplot(wordf,aes(label=type,size=M3,color=color))+
  geom_text_wordcloud_area(shape="circle")+
  scale_size_area(max_size=20)+
  theme_minimal()
ggsave(filename = "regulon_M3.pdf", plot = p1, device = 'pdf', width = 14, height = 10, units = 'cm')
#
wordf <- arrange(Module_mean[,c(1,5)], -M4)
p1=ggplot(wordf,aes(label=type,size=M4,color=color))+
  geom_text_wordcloud_area(shape="circle")+
  scale_size_area(max_size=20)+
  theme_minimal()
ggsave(filename = "regulon_M4.pdf", plot = p1, device = 'pdf', width = 14, height = 10, units = 'cm')
#
wordf <- arrange(Module_mean[,c(1,6)], -M5)
p1=ggplot(wordf,aes(label=type,size=M5,color=color))+
  geom_text_wordcloud_area(shape="circle")+
  scale_size_area(max_size=20)+
  theme_minimal()
ggsave(filename = "regulon_M5.pdf", plot = p1, device = 'pdf', width = 14, height = 10, units = 'cm')
#
wordf <- arrange(Module_mean[,c(1,7)], -M6)
p1=ggplot(wordf,aes(label=type,size=M6,color=color))+
  geom_text_wordcloud_area(shape="circle")+
  scale_size_area(max_size=20)+
  theme_minimal()
ggsave(filename = "regulon_M6.pdf", plot = p1, device = 'pdf', width = 14, height = 10, units = 'cm')

# Figure 2C
# plot the Module-TF-gene network using cytoscape












