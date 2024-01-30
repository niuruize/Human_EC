################################################################################
#---Figure 4 A
################################################################################
# The input data were obtained from cross "EC_integration.R"
load("Human_EC.RData")
table(scRNA$Organ_celltype)
# plot the bar graph using Prism 9

################################################################################
#---Figure 4 B
################################################################################
# library packages
library(Seurat);library(harmony);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl);library(harmony)
# Marker genes for cell subsets in different tissues were calculated separately；
# Take brain for example
load("Brain_Harmony.RData")
table(scRNA$celltype_label2)
scRNA = scRNA[,scRNA$celltype_label2 %in% c("artery","capillary","vein")]
DefaultAssay(scRNA)="RNA"
scRNA <- NormalizeData(scRNA)
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
#鉴定细胞类型的marker gene
Idents(scRNA) <- "celltype_label2"
DefaultAssay(scRNA) <- "RNA"
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = FALSE, test.use = "wilcox")  ##耗时久
write.table(markers,file="markers.txt",quote=F,sep="\t",row.names=F,col.names=T)
# rank and merge
# Take artery for example
adipose <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/adipose/3_conserved _marker/artery.xlsx"))
names(adipose) <- c("gene","adipose")
bladder <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/Bladder/3_conserved _marker/artery.xlsx"))
names(bladder) <- c("gene","bladder")
brain <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/brain/3_conserved _marker/artery.xlsx"))
names(brain) <- c("gene","brain")
esophagus <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/Esophagus/3_conserved _marker/artery.xlsx"))
names(esophagus) <- c("gene","esophagus")
MG <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/gland/3_conserved _marker/artery.xlsx"))
names(MG) <- c("gene","MG")
heart <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/heart/3_conserved _marker/artery.xlsx"))
names(heart) <- c("gene","heart")
intestinal <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/Intestinal/3_conserved _marker/artery.xlsx"))
names(intestinal) <- c("gene","intestinal")
kidney <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/kidney/3_conserved _marker/artery.xlsx"))
names(kidney) <- c("gene","kidney")
liver <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/liver/3_conserved _marker/artery.xlsx"))
names(liver) <- c("gene","liver")
lung <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/lung/3_conserved _marker/artery.xlsx"))
names(lung) <- c("gene","lung")
placenta <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/maternal_fetal/placenta/3_conserved _marker/artery.xlsx"))
names(placenta) <- c("gene","placenta")
muscle <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/muscle/3_conserved _marker/artery.xlsx"))
names(muscle) <- c("gene","muscle")
ovary <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/ovary/3_conserved _marker/artery.xlsx"))
names(ovary) <- c("gene","ovary")
pancreas <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/pancreas/3_conserved _marker/artery.xlsx"))
names(pancreas) <- c("gene","pancreas")
PIT <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/PIT/3_conserved _marker/artery.xlsx"))
names(PIT) <- c("gene","PIT")
SC <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/SC/3_conserved _marker/artery.xlsx"))
names(SC) <- c("gene","SC")
skin <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/skin/3_conserved _marker/artery.xlsx"))
names(skin) <- c("gene","skin")
testis <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/testis/3_conserved _marker/artery.xlsx"))
names(testis) <- c("gene","testis")
thyroid <-  data.frame(readxl::read_xlsx("/Users/niuruize/Downloads/scRNA/EC/thyroid/3_conserved _marker/artery.xlsx"))
names(thyroid) <- c("gene","thyroid")

##
artery <- merge(adipose,bladder, by="gene",all=TRUE)
artery <- merge(artery,brain, by="gene",all=TRUE)
artery <- merge(artery,esophagus, by="gene",all=TRUE)
artery <- merge(artery,heart, by="gene",all=TRUE)
artery <- merge(artery,intestinal, by="gene",all=TRUE)
artery <- merge(artery,kidney, by="gene",all=TRUE)
artery <- merge(artery,liver, by="gene",all=TRUE)
artery <- merge(artery,lung, by="gene",all=TRUE)
artery <- merge(artery,MG, by="gene",all=TRUE)
artery <- merge(artery,muscle, by="gene",all=TRUE)
artery <- merge(artery,ovary, by="gene",all=TRUE)
artery <- merge(artery,pancreas, by="gene",all=TRUE)
artery <- merge(artery,PIT, by="gene",all=TRUE)
artery <- merge(artery,placenta, by="gene",all=TRUE)
artery <- merge(artery,SC, by="gene",all=TRUE)
artery <- merge(artery,skin, by="gene",all=TRUE)
artery <- merge(artery,testis, by="gene",all=TRUE)
artery <- merge(artery,thyroid, by="gene",all=TRUE)

head(artery)
write.csv(artery,"artery_marker.csv", row.names = F)
##
artery <- read.csv('artery_marker.csv',row.names = 1)
# plot the scatter diagram using Prism 9

################################################################################
#---Figure 4 C
################################################################################
library(dplyr)
library(GeneOverlap)
data(GeneOverlap)
#Select variables
variable_list=list("Group")
for (k in 1:length(variable_list)) {
  variable=variable_list[[k]]  
  celltype = c("Adipose_Artery","Adipose_capillary","Adipose_lymphatic","Adipose_vein",
               "Bladder_Artery","Bladder_capillary","Bladder_lymphatic","Bladder_vein",
               "Brain_Artery","Brain_capillary","Brain_vein",
               "Esophagus_artery","Esophagus_capillary","Esophagus_vein",
               "MG_artery","MG_capillary","MG_vein",
               "Heart_Artery","Heart_capillary","Heart_lymphatic","Heart_vein",
               "Intestinal_Artery","Intestinal_lymphatic","Intestinal_vein",
               "Kidney_Artery","Kidney_capillary","Kidney_vein",
               "Liver_Artery","Liver_capillary","Liver_lymphatic","Liver_vein",
               "Lung_Artery","Lung_capillary","Lung_vein",
               "Placenta_Artery","Placenta_capillary","Placenta_vein",
               "Muscle_Artery","Muscle_capillary","Muscle_lymphatic","Muscle_vein",
               "Ovary_Artery","Ovary_capillary","Ovary_lymphatic","Ovary_vein",
               "Pancreas_Artery","Pancreas_capillary","Pancreas_lymphatic","Pancreas_vein",
               "PIT_Artery","PIT_vein",
               "SC_Artery","SC_capillary","SC_vein",
               "Skin_Artery","Skin_lymphatic","Skin_vein",
               "Testis_Artery","Testis_lymphatic","Testis_vein",
               "Thyroid_Artery","Thyroid_capillary","Thyroid_lymphatic","Thyroid_vein",
               "Uterus_capillary","Uterus_vein")
  #List significantly upregulated genes
  
  Marker_genes_list<-list()
  for (i in 1:length(celltype)) {
    ct1=celltype[i]
    ct1=gsub("/","_",ct1)
    input_file=paste0("/Users/niuruize/Downloads/scRNA/EC/EC/4_marker/all/",ct1,".xlsx")
    topgenes=readxl::read_xlsx(input_file)
    Marker_genes_full_significant=topgenes
    #Marker_genes_full_significant=subset(Marker_genes_full,fdr < 0.05)
    #Marker_genes_full_significant=subset(Marker_genes_full_significant,logFC > 0)
    Marker_genes=Marker_genes_full_significant[,1:2]
    Marker_genes[] <- lapply(Marker_genes, as.character)
    if (nrow(Marker_genes)>0){
      Marker_genes_list[[i]]<-Marker_genes[[1]]
      names(Marker_genes_list)[i]<- celltype[i]
    }
  }
  
  #Remove empty elements
  Marker_genes_list=Marker_genes_list[lengths(Marker_genes_list) != 0]
  
  for (c in 1:length(celltype)) {
    ct1=celltype[[c]]
    ct1_name=gsub("/","_",ct1)
    for (b in 1:length(celltype)){
      ct2=celltype[[b]]
      ct2_name=gsub("/","_",ct2)
      
      input_file1=paste0("/Users/niuruize/Downloads/scRNA/EC/EC/4_marker/all/",ct1_name,"_Aging_Adult_DEGs",".csv")
      input_file2=paste0("/Users/niuruize/Downloads/scRNA/EC/EC/4_marker/all/",ct2_name,"_Aging_Adult_DEGs",".csv")
      
      topgenes1=read.csv(input_file1,header=TRUE,row.names=1)
      
      topgenes2=read.csv(input_file2,header=TRUE,row.names=1)
      
      All_genes=full_join(topgenes1,topgenes2,by="symbol")
      
      genome_size=nrow(All_genes)
      
      gom.obj <- newGOM(Marker_genes_list,Marker_genes_list, genome.size = genome_size)
      GeneOverlap_Matrix=getMatrix(gom.obj, name="pval")
      GeneOverlap_Matrix_df=as.data.frame(GeneOverlap_Matrix)
      
      myData=GeneOverlap_Matrix_df[ct1,ct2,drop=FALSE]
      
      if (b==1){
        results_int=myData
      }
      if(b>1){
        results_int=cbind(results_int,myData)
      }
    }
    if (c==1){
      results_table=results_int
    }
    if (c>1){
      results_table=rbind(results_table,results_int)
    }
  }
  
  n_rows=nrow(results_table)
  n_columns=ncol(results_table)
  factor=n_rows*n_columns
  pvalue_adjusted_matrix <- apply(results_table,2,p.adjust,method="bonferroni", n = factor)
  pvalue_adjusted_matrix_minusLog10=-log10(pvalue_adjusted_matrix)
  pvalue_adjusted_matrix_minusLog10[is.infinite(pvalue_adjusted_matrix_minusLog10)] <- 300
  
  setwd("/Users/niuruize/Downloads/scRNA/EC/EC/4_marker/all")
  output_file=paste0("EC_marker_P",".csv")
  write.csv(pvalue_adjusted_matrix_minusLog10,file=output_file)
  
  
  #GeneOverlap_Matrix=getMatrix(gom.obj, name="odds.ratio")
  GeneOverlap_Matrix=getMatrix(gom.obj, name="Jaccard")
  #GeneOverlap_Matrix[is.infinite(GeneOverlap_Matrix)] <- 20
  odds_ratio_Matrix = as.data.frame(GeneOverlap_Matrix)
  output_file=paste0("EC_marker_OR",".csv")
  write.csv(odds_ratio_Matrix,file=output_file)
}

#
# https://www.sohu.com/a/713708194_121123706 
library(pheatmap)
#r2 <- odds_ratio_Matrix
#p2 <- pvalue_adjusted_matrix_minusLog10

r2 <- read.csv('EC_marker_OR.csv',header=TRUE,row.names=1)
p2 <- read.csv('EC_marker_P.csv',header=TRUE,row.names=1)

r2[r2 >= 10] <- 10

#
p2[p2 <= 10] <- ""
p2[p2 >= 100] <- "**"
p2[p2 > 10] <- "*"

#mycol<-colorRampPalette(c( "#0f86a9", "white", "#ed8b10"))(200)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
mycol<-colorRampPalette(c("white", "#0f86a9"))(200)
#
p1 = pheatmap(r2,scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=12,fontsize_row=8,fontsize_col=8,cellwidth=15,
              cellheight=15,cluster_rows=T,cluster_cols=T,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "EC_marker_overlap3.pdf", p1,device = 'pdf', width = 40, height = 40, units = 'cm')
p1 = pheatmap(r2,scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=12,fontsize_row=8,fontsize_col=8,cellwidth=15,
              cellheight=15,cluster_rows=F,cluster_cols=F,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "EC_aging_overlap_up4.pdf", p1,device = 'pdf', width = 40, height = 40, units = 'cm')










