################################################################################
#---Figure 5 scFEA
################################################################################
#loading R packages
library(Seurat);library(ggplot2)
Idents(scRNA) <- "Organ"
table(scRNA$Organ)
scRNA$Organ <- as.factor(as.character(scRNA$Organ))
scRNA1 <- subset(scRNA, subset = nFeature_RNA > 500)
table(scRNA1$Organ)
save(scRNA1,file="EC_500_all.RData")
# Extract expression matrix
EC20000 <- as.matrix(scRNA1@assays$RNA@data[,1:20000])
write.csv(EC20000, file='EC20000.csv', row.names = T)
EC40000 <- as.matrix(scRNA1@assays$RNA@data[,20001:40000])
write.csv(EC40000, file='EC40000.csv', row.names = T)
EC60000 <- as.matrix(scRNA1@assays$RNA@data[,40001:60000])
write.csv(EC60000, file='EC60000.csv', row.names = T)
EC80000 <- as.matrix(scRNA1@assays$RNA@data[,60001:80000])
write.csv(EC80000, file='EC80000.csv', row.names = T)
EC100000 <- as.matrix(scRNA1@assays$RNA@data[,80001:100000])
write.csv(EC100000, file='EC100000.csv', row.names = T)
EC120000 <- as.matrix(scRNA1@assays$RNA@data[,100001:120000])
write.csv(EC120000, file='EC120000.csv', row.names = T)
EC140000 <- as.matrix(scRNA1@assays$RNA@data[,120001:140000])
write.csv(EC140000, file='EC140000.csv', row.names = T)
EC160000 <- as.matrix(scRNA1@assays$RNA@data[,140001:160000])
write.csv(EC160000, file='EC160000.csv', row.names = T)
EC180000 <- as.matrix(scRNA1@assays$RNA@data[,160001:188037])
write.csv(EC180000, file='EC180000.csv', row.names = T)
# scFEA analysis used python
nohup python src/scFEA.py --data_dir data --input_dir input \
                    --test_file EC20000.csv \
                    --moduleGene_file module_gene_m168.csv \
                    --stoichiometry_matrix cmMat_c70_m168.csv \
                    --output_flux_file output/EC20000_flux.csv \
                    --output_balance_file output/EC20000_balance.csv &
################################################################################
#---Figure 5A
################################################################################
####
predFlux1 <- read.csv('output/EC20000_flux.csv', header = T, row.names = 1)
predFlux2 <- read.csv('output/EC40000_flux.csv', header = T, row.names = 1)
predFlux3 <- read.csv('output/EC60000_flux.csv', header = T, row.names = 1)
predFlux4 <- read.csv('output/EC80000_flux.csv', header = T, row.names = 1)
predFlux5 <- read.csv('output/EC100000_flux.csv', header = T, row.names = 1)
predFlux6 <- read.csv('output/EC120000_flux.csv', header = T, row.names = 1)
predFlux7 <- read.csv('output/EC140000_flux.csv', header = T, row.names = 1)
predFlux8 <- read.csv('output/EC160000_flux.csv', header = T, row.names = 1)
predFlux9 <- read.csv('output/EC180000_flux.csv', header = T, row.names = 1)

predFlux <- rbind(predFlux1,predFlux2,predFlux3,predFlux4,predFlux5,predFlux6,predFlux7,predFlux8,predFlux9)
scRNA <- AddMetaData(obj, metadata = predFlux)
Organ_scFEA<-FetchData(scRNA, vars = c("Organ",colnames(predFlux)))
Organ_scFEA_mean <- aggregate(Organ_scFEA[,2:169], by=list(type=Organ_scFEA$Organ),mean)
rownames(Organ_scFEA_mean) <- Organ_scFEA_mean$type
Organ_scFEA_mean <- Organ_scFEA_mean[,-1]
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
Human_M168 <- read.csv('Human_M168_information.symbols.csv', header = T, row.names = 1)
#annotation_row <- data.frame(Metabolites = Human_M168$name)
#rownames(annotation_row) <- rownames(Human_M168)
p1 = pheatmap(t(Organ_scFEA_mean),show_colnames = T,scale = "row", cluster_cols = T,
              treeheight_col = 10,color= mycol, treeheight_row=15,annotation_row = annotation_row)
ggsave(filename = "scFEA_flux500_all_Organ2.pdf", p1,device = 'pdf', width = 20, height = 60, units = 'cm')

################################################################################
#---Figure 5B
################################################################################
####
predFlux1 <- read.csv('output/EC20000_balance.csv', header = T, row.names = 1)
predFlux2 <- read.csv('output/EC40000_balance.csv', header = T, row.names = 1)
predFlux3 <- read.csv('output/EC60000_balance.csv', header = T, row.names = 1)
predFlux4 <- read.csv('output/EC80000_balance.csv', header = T, row.names = 1)
predFlux5 <- read.csv('output/EC100000_balance.csv', header = T, row.names = 1)
predFlux6 <- read.csv('output/EC120000_balance.csv', header = T, row.names = 1)
predFlux7 <- read.csv('output/EC140000_balance.csv', header = T, row.names = 1)
predFlux8 <- read.csv('output/EC160000_balance.csv', header = T, row.names = 1)
predFlux9 <- read.csv('output/EC180000_balance.csv', header = T, row.names = 1)

predFlux <- rbind(predFlux1,predFlux2,predFlux3,predFlux4,predFlux5,predFlux6,predFlux7,predFlux8,predFlux9)
scRNA <- AddMetaData(scRNA1, metadata = predFlux)
Organ_scFEA<-FetchData(scRNA, vars = c("Organ",colnames(predFlux)))
Organ_scFEA_mean <- aggregate(Organ_scFEA[,2:71], by=list(type=Organ_scFEA$Organ),mean)
rownames(Organ_scFEA_mean) <- Organ_scFEA_mean$type
Organ_scFEA_mean <- Organ_scFEA_mean[,-1]
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap(t(Organ_scFEA_mean),show_colnames = T,scale = "row", cluster_cols = T,treeheight_col = 10, treeheight_row=15,color= mycol)
ggsave(filename = "scFEA_balance_500_all_Organ.pdf", p1,device = 'pdf', width = 18, height = 30, units = 'cm')

################################################################################
#---Figure 5C
################################################################################
markers <- read_xlsx("SLC_gene.xlsx")
p1 = jjDotPlot(scRNA, gene = SLC_gene$gene,id = 'Organ',xtree = T,ytree = T,rescale = T,rescale.min = 0,rescale.max = 1)
ggsave(filename = "EC_SLC_marker_heatmap.pdf", plot = p1, device = 'pdf', width = 21, height = 21, units = 'cm')


