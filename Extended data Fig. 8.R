################################################################################
# Extended data Fig. 8A
################################################################################
# Extract expression matrix
library(Seurat);library(ggplot2)
load('MacFas_EC.RData')
Idents(scRNA) <- "Organ"
table(scRNA$Organ)
scRNA$Organ <- as.factor(as.character(scRNA$Organ))
scRNA1 <- subset(scRNA, subset = nFeature_RNA > 500)
table(scRNA1$Organ)
save(scRNA1,file="EC_500_all.RData")
macaque_EC <- as.matrix(scRNA1@assays$RNA@data)
write.csv(macaque_EC, file='macaque_EC.csv', row.names = T)
# scFEA
nohup python src/scFEA.py --data_dir data --input_dir input \
                    --test_file macaque_EC.csv \
                    --moduleGene_file module_gene_m168.csv \
                    --stoichiometry_matrix cmMat_c70_m168.csv \
                    --output_flux_file output/macaque_EC_flux.csv \
                    --output_balance_file output/macaque_EC_balance.csv &
# visualization
predFlux <- read.csv('output/macaque_EC_flux.csv', header = T, row.names = 1)
scRNA <- AddMetaData(obj, metadata = predFlux)
Human_M168 <- read.csv('Human_M168_information.symbols.csv', header = T, row.names = 1)
annotation_row <- data.frame(Metabolites = Human_M168$name)
rownames(annotation_row) <- rownames(Human_M168)
Organ_scFEA<-FetchData(scRNA, vars = c("Organ",colnames(predFlux)))
Organ_scFEA_mean <- aggregate(Organ_scFEA[,2:169], by=list(type=Organ_scFEA$Organ),mean)
rownames(Organ_scFEA_mean) <- Organ_scFEA_mean$type
Organ_scFEA_mean <- Organ_scFEA_mean[,-1]
Organ_scFEA_mean = Organ_scFEA_mean[apply(Organ_scFEA_mean, 1, function(x) sd(x)!=0),] 
Organ_scFEA_mean = Organ_scFEA_mean[,apply(Organ_scFEA_mean, 2, function(x) sd(x)!=0)]
p1 = pheatmap(t(Organ_scFEA_mean),show_colnames = T,scale = "row", cluster_cols = T,
              treeheight_col = 10,color= mycol, treeheight_row=15,annotation_row = annotation_row)
ggsave(filename = "scFEA_flux500_all_Organ1.pdf", p1,device = 'pdf', width = 16, height = 60, units = 'cm')
colnames(Organ_scFEA_mean)=Human_M168$ID_name

################################################################################
# Extended data Fig. 8B
################################################################################
# Extract expression matrix
library(Seurat);library(ggplot2)
load('Mouse_EC.RData')
Idents(scRNA) <- "Organ"
table(scRNA$Organ)
scRNA$Organ <- as.factor(as.character(scRNA$Organ))
scRNA1 <- subset(scRNA, subset = nFeature_RNA > 500)
table(scRNA1$Organ)
save(scRNA1,file="EC_500_all.RData")
Mouse_EC <- as.matrix(scRNA1@assays$RNA@data)
write.csv(Mouse_EC, file='Mouse_EC.csv', row.names = T)
# scFEA
nohup python src/scFEA.py --data_dir data --input_dir input \
                    --test_file Mouse_EC.csv \
                    --moduleGene_file module_gene_complete_mouse_m168.csv \
                    --stoichiometry_matrix cmMat_c70_m168.csv \
                    --output_flux_file output/Mouse_EC_flux.csv \
                    --output_balance_file output/Mouse_EC_balance.csv &
# visualization
predFlux <- read.csv('output/Mouse_EC_balance.csv', header = T, row.names = 1)
scRNA <- AddMetaData(obj, metadata = predFlux)
Human_M168 <- read.csv('Human_M168_information.symbols.csv', header = T, row.names = 1)
annotation_row <- data.frame(Metabolites = Human_M168$name)
rownames(annotation_row) <- rownames(Human_M168)
Organ_scFEA<-FetchData(scRNA, vars = c("Organ",colnames(predFlux)))
Organ_scFEA_mean <- aggregate(Organ_scFEA[,2:169], by=list(type=Organ_scFEA$Organ),mean)
rownames(Organ_scFEA_mean) <- Organ_scFEA_mean$type
Organ_scFEA_mean <- Organ_scFEA_mean[,-1]
Organ_scFEA_mean = Organ_scFEA_mean[apply(Organ_scFEA_mean, 1, function(x) sd(x)!=0),] 
Organ_scFEA_mean = Organ_scFEA_mean[,apply(Organ_scFEA_mean, 2, function(x) sd(x)!=0)]
p1 = pheatmap(t(Organ_scFEA_mean),show_colnames = T,scale = "row", cluster_cols = T,
              treeheight_col = 10,color= mycol, treeheight_row=15,annotation_row = annotation_row)
ggsave(filename = "scFEA_flux500_all_Organ1.pdf", p1,device = 'pdf', width = 15, height = 60, units = 'cm')
colnames(Organ_scFEA_mean)=Human_M168$ID_name





