
################################################################################
# Extended data Fig. 10A
################################################################################
#scTenifoldKnk
library(scTenifoldKnk);library(Seurat)
load('Liver_Endo')
table(scRNA$Group)
scRNA_seq <- scRNA[,scRNA$Group %in% c("Adult")]
scRNA_seq$Group <- as.factor(as.character(scRNA_seq$Group))
X <- data.frame(scRNA_seq@assays$RNA@counts)
save(X ,file="Liver_Adult_KO_counts.RData")
##
library(scTenifoldKnk);library(Seurat);library(Matrix)
args = commandArgs(trailingOnly=TRUE)
scQC <- function(X, mtThreshold = 0.1, minLSize = 1000){
  require(Matrix)
  if(class(X) == 'Seurat'){
    countMatrix <- X@assays$RNA@counts
  } else {
    countMatrix <- X
  }
  librarySize <- colSums(countMatrix)
  countMatrix <- countMatrix[,librarySize >= minLSize]
  librarySize <- colSums(countMatrix)
  mtGenes <- grep('^MT-',toupper(rownames(countMatrix)))
  nGenes <- colSums(countMatrix != 0)
  
  genesLM <- lm(nGenes~librarySize)
  genesLM <- as.data.frame(predict(genesLM, data.frame(librarySize), interval = 'prediction'))
  
  if(isTRUE(length(mtGenes) > 0)){
    mtCounts <- colSums(countMatrix[grep('^MT-',toupper(rownames(countMatrix))),])
    mtProportion <- mtCounts/librarySize
    mtLM <- lm(mtCounts~librarySize)
    mtLM <- as.data.frame(predict(mtLM, data.frame(librarySize), interval = 'prediction'))
    selectedCells <- ((mtCounts > mtLM$lwr) & (mtCounts < mtLM$upr) & (nGenes > genesLM$lwr) & (nGenes < genesLM$upr) & (mtProportion <= mtThreshold) & (librarySize < 2 * mean(librarySize)))
  } else {
    selectedCells <- ((nGenes > genesLM$lwr) & (nGenes < genesLM$upr) & (librarySize < 2 * mean(librarySize)))
  }
  selectedCells <- colnames(countMatrix)[selectedCells]
  if(class(X) == 'Seurat'){
    X <- subset(X, cells = selectedCells)
  } else {
    X <- countMatrix[,selectedCells]
  }
  return(X)
}

# Running scTenifoldKnk
load("Liver_Adult_KO_counts.RData")
X <- scQC(X)
tempX <- X[,sample(seq_len(ncol(X)),6000)]
tempX <- tempX[rowMeans(X != 0) > 0.05,]
oX <- scTenifoldKnk(countMatrix = tempX, gKO = 'NUFIP2',qc = TRUE, nCores = 16)
save(oX ,file="Liver_NUFIP2_counts.RData")
write.csv(oX$diffRegulation, file = paste0('NUFIP2_counts_Stability_',args[1], '.csv'), row.names = FALSE)
##
library(fgsea);library(ggplot2);library(enrichR);library(igraph)
source('/Users/niuruize/Downloads/scRNA/R_code/plotDR.R')
source('/Users/niuruize/Downloads/scRNA/R_code/plotKO.R')
source('/Users/niuruize/Downloads/scRNA/R_code/hsa2mmu_SYMBOL.R')
load('Liver_NUFIP2_counts.RData')
Liver_NUFIP2 = oX
write.csv(Liver_NUFIP2$diffRegulation, 'Liver_NUFIP2_counts.csv')
dGenes <- Liver_NUFIP2$diffRegulation$gene[Liver_NUFIP2$diffRegulation$p.adj < 0.05]

#png('dr2_Liver_NUFIP2.png', width = 2000, height = 2000, res = 300)
pdf("dr2_Liver_NUFIP2.pdf", width = 5,height = 5)
plotDR(Liver_NUFIP2)
dev.off()

#png('ego2_Liver_NUFIP2.png', width = 3000, height = 3000, res = 300, bg = NA)
pdf("ego2_Liver_NUFIP2.pdf", width = 10,height = 10)
X <- Liver_NUFIP2
gKO <- 'NUFIP2'
q <- 0.995
gList <- unique(c(gKO, X$diffRegulation$gene[X$diffRegulation$p.adj < 0.05]))
sCluster <- as.matrix(X$tensorNetworks$WT[gList,gList])
koInfo <- sCluster[gKO,]
sCluster[abs(sCluster) <= quantile(abs(sCluster), q)] <- 0
sCluster[gKO,] <- koInfo
diag(sCluster) <- 0
sCluster <-  reshape2::melt(as.matrix(sCluster))
colnames(sCluster) <- c('from', 'to', 'W')
sCluster <- sCluster[sCluster$W != 0,]
netPlot <- graph_from_data_frame(sCluster, directed = TRUE)
dPlot <- centr_degree(netPlot)$res
W <- rep(1,nrow(sCluster))
sG   <- (names(V(netPlot))[dPlot > 1])[-1]
W[sCluster$from %in% sG] <- 0.2
W[sCluster$to %in% sG] <- 0.2
W[sCluster$from %in% gKO] <- 1
W[sCluster$from %in% gKO & sCluster$to %in% sG] <- 0.8
set.seed(1)
layPlot <- layout_with_fr(netPlot, weights = W)
dPlot <- (dPlot/max(dPlot))*20
dbs <- listEnrichrDbs()
E <- enrichr(gList, c("GO_Biological_Process_2023","KEGG_2021_Human"))
E <- do.call(rbind.data.frame, E)
#E <- E[E$Adjusted.P.value < 0.05,]
E <- E[order(E$Adjusted.P.value),]
E$Term <- unlist(lapply(strsplit(E$Term,''), function(X){
  X[1] <- toupper(X[1])
  X <- paste0(X,collapse = '')
  X <- gsub('\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+','',X)
  X <- gsub("'s",'',X)
  X <- unlist(strsplit(X,','))[1]
  X <- gsub('[[:blank:]]$','',X)
  return(X)
}))
write.csv(E,"NUFIP2_EGO.csv") 
E <- E[E$Term %in% c('Inflammatory Response','Positive Regulation Of Phagocytosis','Complement Activation','Regulation Of Vasoconstriction','Glycolipid Transport'),]
#E <- E[c(1,2,5,7),]
tPlot <- strsplit(E$Genes, ';')
pPlot <- matrix(0,nrow = length(V(netPlot)), ncol = nrow(E))
rownames(pPlot) <- toupper(names(V(netPlot)))
for(i in seq_along(tPlot)){
  pPlot[unlist(tPlot[i]),i] <- 1
}
pPlot <- lapply(seq_len(nrow(pPlot)), function(X){as.vector(pPlot[X,])})
names(pPlot) <- names(V(netPlot))
tPlot <- unique(unlist(tPlot))
eGenes <- toupper(names(V(netPlot))) %in% tPlot
vColor <- rgb(195/255, 199/255, 198/255 ,0.3)
pieColors <- list(hcl.colors(nrow(E), palette = 'Zissou 1', alpha = 0.7))
par(mar=c(4,0,0,0), xpd = TRUE)
suppressWarnings(plot(netPlot,
                      layout = layPlot,
                      edge.arrow.size=.2,
                      vertex.label.color="black",
                      vertex.shape = ifelse(eGenes,'pie','circle'),
                      vertex.pie = pPlot,
                      vertex.size = 10+dPlot,
                      vertex.pie.color=pieColors,
                      vertex.label.family="Times",
                      vertex.label.font=ifelse(eGenes,2,1),
                      edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue'),
                      edge.curved = ifelse(W == 0.2, 0, 0.1),
                      vertex.color = vColor,
                      vertex.frame.color = NA))
sigLevel <- formatC(E$Adjusted.P.value, digits = 2, format = 'g', width = 0, drop0trailing = TRUE)
gSetNames <- lengths(strsplit(E$Genes, ';'))
gSetNames <- paste0('(', gSetNames,') ', E$Term, ' FDR = ', sigLevel)
legend(x = -1.05, y = -1.05, legend = gSetNames, bty = 'n', ncol = 2, cex = 1, col = unlist(pieColors), pch = 16)
dev.off()


zLiver_NUFIP2 <- Liver_NUFIP2$diffRegulation$Z
names(zLiver_NUFIP2) <- toupper(Liver_NUFIP2$diffRegulation$gene)

MGI <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human')
set.seed(1)
E <- fgseaMultilevel(MGI, zLiver_NUFIP2)

#png('gsea1_Liver_NUFIP2.png', width = 1000, height = 1000, res = 300)
pdf("gsea1_Liver_NUFIP2.pdf", width = 5,height = 5)
gSet <- 'Alzheimer disease'
plotEnrichment(MGI[[gSet]], zLiver_NUFIP2) +
  labs(
    title = 'Alzheimer disease',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea2_Liver_NUFIP2.png', width = 1000, height = 1000, res = 300)
pdf("gsea2_Liver_NUFIP2.pdf", width = 5,height = 5)
gSet <- 'Non-alcoholic fatty liver disease (NAFLD)'
plotEnrichment(MGI[[gSet]], zLiver_NUFIP2) +
  labs(
    title = 'Non-alcoholic fatty liver disease (NAFLD)',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea3_Liver_NUFIP2.png', width = 1000, height = 1000, res = 300)
pdf("gsea3_Liver_NUFIP2.pdf", width = 5,height = 5)
gSet <- 'Parkinson disease'
plotEnrichment(MGI[[gSet]], zLiver_NUFIP2) +
  labs(
    title = 'Parkinson disease',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score')  + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea3_Liver_NUFIP2.png', width = 1000, height = 1000, res = 300)
pdf("gsea4_Liver_NUFIP2.pdf", width = 5,height = 5)
gSet <- 'Huntington disease'
plotEnrichment(MGI[[gSet]], zLiver_NUFIP2) +
  labs(
    title = 'Huntington disease',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score')  + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea3_Liver_NUFIP2.png', width = 1000, height = 1000, res = 300)
pdf("gsea5_Liver_NUFIP2.pdf", width = 5,height = 5)
gSet <- 'Fatty acid degradation'
plotEnrichment(MGI[[gSet]], zLiver_NUFIP2) +
  labs(
    title = 'Fatty acid degradation',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score')  + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea3_Liver_NUFIP2.png', width = 1000, height = 1000, res = 300)
pdf("gsea6_Liver_NUFIP2.pdf", width = 5,height = 5)
gSet <- 'Arginine and proline metabolism'
plotEnrichment(MGI[[gSet]], zLiver_NUFIP2) +
  labs(
    title = 'Arginine and proline metabolism',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score')  + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea3_Liver_NUFIP2.png', width = 1000, height = 1000, res = 300)
pdf("gsea7_Liver_NUFIP2.pdf", width = 5,height = 5)
gSet <- 'DNA replication'
plotEnrichment(MGI[[gSet]], zLiver_NUFIP2) +
  labs(
    title = 'DNA replication',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score')  + theme(plot.title = element_text(face = 2, size = 25))
dev.off()
################################################################################
# Extended data Fig. 10B
################################################################################
#scTenifoldKnk
load(Bladder_Endo)
table(scRNA$Group)
scRNA_seq <- scRNA[,scRNA$Group %in% c("Aging")]
scRNA_seq$Group <- as.factor(as.character(scRNA_seq$Group))
X <- data.frame(scRNA_seq@assays$RNA@counts)
save(X ,file="Bladder_Aging_KO_counts.RData")
##
library(scTenifoldKnk);library(Seurat);library(Matrix)
args = commandArgs(trailingOnly=TRUE)
scQC <- function(X, mtThreshold = 0.1, minLSize = 1000){
  require(Matrix)
  if(class(X) == 'Seurat'){
    countMatrix <- X@assays$RNA@counts
  } else {
    countMatrix <- X
  }
  librarySize <- colSums(countMatrix)
  countMatrix <- countMatrix[,librarySize >= minLSize]
  librarySize <- colSums(countMatrix)
  mtGenes <- grep('^MT-',toupper(rownames(countMatrix)))
  nGenes <- colSums(countMatrix != 0)
  
  genesLM <- lm(nGenes~librarySize)
  genesLM <- as.data.frame(predict(genesLM, data.frame(librarySize), interval = 'prediction'))
  
  if(isTRUE(length(mtGenes) > 0)){
    mtCounts <- colSums(countMatrix[grep('^MT-',toupper(rownames(countMatrix))),])
    mtProportion <- mtCounts/librarySize
    mtLM <- lm(mtCounts~librarySize)
    mtLM <- as.data.frame(predict(mtLM, data.frame(librarySize), interval = 'prediction'))
    selectedCells <- ((mtCounts > mtLM$lwr) & (mtCounts < mtLM$upr) & (nGenes > genesLM$lwr) & (nGenes < genesLM$upr) & (mtProportion <= mtThreshold) & (librarySize < 2 * mean(librarySize)))
  } else {
    selectedCells <- ((nGenes > genesLM$lwr) & (nGenes < genesLM$upr) & (librarySize < 2 * mean(librarySize)))
  }
  selectedCells <- colnames(countMatrix)[selectedCells]
  if(class(X) == 'Seurat'){
    X <- subset(X, cells = selectedCells)
  } else {
    X <- countMatrix[,selectedCells]
  }
  return(X)
}

# Running scTenifoldKnk
load("Bladder_Aging_KO_counts.RData")
X <- scQC(X)
tempX <- X[,sample(seq_len(ncol(X)),300)]
tempX <- tempX[rowMeans(X != 0) > 0.05,]
oX <- scTenifoldKnk(countMatrix = tempX, gKO = 'TMSB4X',qc = FALSE, nCores = 16)
save(oX ,file="Bladder_TMSB4X_data.RData")
write.csv(oX$diffRegulation, file = paste0('TMSB4X_data_Stability_',args[1], '.csv'), row.names = FALSE)
##
library(fgsea);library(ggplot2);library(enrichR);library(igraph)
source('/Users/niuruize/Downloads/scRNA/R_code/plotDR.R')
source('/Users/niuruize/Downloads/scRNA/R_code/plotKO.R')
source('/Users/niuruize/Downloads/scRNA/R_code/hsa2mmu_SYMBOL.R')
load('Bladder_TMSB4X_data.RData')
Bladder_TMSB4X = oX
write.csv(Bladder_TMSB4X$diffRegulation, 'Bladder_TMSB4X_counts.csv')
dGenes <- Bladder_TMSB4X$diffRegulation$gene[Bladder_TMSB4X$diffRegulation$p.adj < 0.05]
#png('dr2_Bladder_TMSB4X.png', width = 2000, height = 2000, res = 300)
pdf("dr2_Bladder_TMSB4X.pdf", width = 5,height = 5)
plotDR(Bladder_TMSB4X)
dev.off()

#png('ego2_Bladder_TMSB4X.png', width = 3000, height = 3000, res = 300, bg = NA)
pdf("ego2_Bladder_TMSB4X.pdf", width = 10,height = 10)
X <- Bladder_TMSB4X
gKO <- 'TMSB4X'
q <- 0.995
gList <- unique(c(gKO, X$diffRegulation$gene[X$diffRegulation$p.adj < 0.05]))
sCluster <- as.matrix(X$tensorNetworks$WT[gList,gList])
koInfo <- sCluster[gKO,]
sCluster[abs(sCluster) <= quantile(abs(sCluster), q)] <- 0
sCluster[gKO,] <- koInfo
diag(sCluster) <- 0
sCluster <-  reshape2::melt(as.matrix(sCluster))
colnames(sCluster) <- c('from', 'to', 'W')
sCluster <- sCluster[sCluster$W != 0,]
netPlot <- graph_from_data_frame(sCluster, directed = TRUE)
dPlot <- centr_degree(netPlot)$res
W <- rep(1,nrow(sCluster))
sG   <- (names(V(netPlot))[dPlot > 1])[-1]
W[sCluster$from %in% sG] <- 0.2
W[sCluster$to %in% sG] <- 0.2
W[sCluster$from %in% gKO] <- 1
W[sCluster$from %in% gKO & sCluster$to %in% sG] <- 0.8
set.seed(1)
layPlot <- layout_with_fr(netPlot, weights = W)
dPlot <- (dPlot/max(dPlot))*20
dbs <- listEnrichrDbs()
E <- enrichr(gList, c("GO_Biological_Process_2023","KEGG_2021_Human"))
E <- do.call(rbind.data.frame, E)
#E <- E[E$Adjusted.P.value < 0.05,]
E <- E[order(E$Adjusted.P.value),]
E$Term <- unlist(lapply(strsplit(E$Term,''), function(X){
  X[1] <- toupper(X[1])
  X <- paste0(X,collapse = '')
  X <- gsub('\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+','',X)
  X <- gsub("'s",'',X)
  X <- unlist(strsplit(X,','))[1]
  X <- gsub('[[:blank:]]$','',X)
  return(X)
}))
write.csv(E,"TMSB4X_EGO.csv")
E <- E[E$Term %in% c('Inflammatory Response','Lipid and atherosclerosis','Regulation Of Vascular Wound Healing','Regulation Of Endothelial Cell Apoptotic Process'),]
#E <- E[c(1,2,5,7),]
tPlot <- strsplit(E$Genes, ';')
pPlot <- matrix(0,nrow = length(V(netPlot)), ncol = nrow(E))
rownames(pPlot) <- toupper(names(V(netPlot)))
for(i in seq_along(tPlot)){
  pPlot[unlist(tPlot[i]),i] <- 1
}
pPlot <- lapply(seq_len(nrow(pPlot)), function(X){as.vector(pPlot[X,])})
names(pPlot) <- names(V(netPlot))
tPlot <- unique(unlist(tPlot))
eGenes <- toupper(names(V(netPlot))) %in% tPlot
vColor <- rgb(195/255, 199/255, 198/255 ,0.3)
pieColors <- list(hcl.colors(nrow(E), palette = 'Zissou 1', alpha = 0.7))
par(mar=c(4,0,0,0), xpd = TRUE)
suppressWarnings(plot(netPlot,
                      layout = layPlot,
                      edge.arrow.size=.2,
                      vertex.label.color="black",
                      vertex.shape = ifelse(eGenes,'pie','circle'),
                      vertex.pie = pPlot,
                      vertex.size = 10+dPlot,
                      vertex.pie.color=pieColors,
                      vertex.label.family="Times",
                      vertex.label.font=ifelse(eGenes,2,1),
                      edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue'),
                      edge.curved = ifelse(W == 0.2, 0, 0.1),
                      vertex.color = vColor,
                      vertex.frame.color = NA))
sigLevel <- formatC(E$Adjusted.P.value, digits = 2, format = 'g', width = 0, drop0trailing = TRUE)
gSetNames <- lengths(strsplit(E$Genes, ';'))
gSetNames <- paste0('(', gSetNames,') ', E$Term, ' FDR = ', sigLevel)
legend(x = -1.05, y = -1.05, legend = gSetNames, bty = 'n', ncol = 2, cex = 1, col = unlist(pieColors), pch = 16)
dev.off()


zBladder_TMSB4X <- Bladder_TMSB4X$diffRegulation$Z
names(zBladder_TMSB4X) <- toupper(Bladder_TMSB4X$diffRegulation$gene)

MGI <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human')
set.seed(1)
E <- fgseaMultilevel(MGI, zBladder_TMSB4X)

#png('gsea1_Bladder_TMSB4X.png', width = 1000, height = 1000, res = 300)
pdf("gsea1_Bladder_TMSB4X.pdf", width = 5,height = 5)
gSet <- 'IL-17 signaling pathway'
plotEnrichment(MGI[[gSet]], zBladder_TMSB4X) +
  labs(
    title = 'IL-17 signaling pathway',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea2_Bladder_TMSB4X.png', width = 1000, height = 1000, res = 300)
pdf("gsea2_Bladder_TMSB4X.pdf", width = 5,height = 5)
gSet <- 'TNF signaling pathway'
plotEnrichment(MGI[[gSet]], zBladder_TMSB4X) +
  labs(
    title = 'TNF signaling pathway',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea3_Bladder_TMSB4X.png', width = 1000, height = 1000, res = 300)
pdf("gsea3_Bladder_TMSB4X.pdf", width = 5,height = 5)
gSet <- 'Apoptosis'
plotEnrichment(MGI[[gSet]], zBladder_TMSB4X) +
  labs(
    title = 'Apoptosis',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score')  + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea3_Bladder_TMSB4X.png', width = 1000, height = 1000, res = 300)
pdf("gsea4_Bladder_TMSB4X.pdf", width = 5,height = 5)
gSet <- 'Necroptosis'
plotEnrichment(MGI[[gSet]], zBladder_TMSB4X) +
  labs(
    title = 'Necroptosis',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score')  + theme(plot.title = element_text(face = 2, size = 25))
dev.off()
write.csv(data.frame(lapply(E, as.character)),"TMSB4X_GSEA.csv") #保存结果

################################################################################
# Extended data Fig. 10C
################################################################################
library(Seurat);library(tidyverse);library(cowplot);library(Matrix);library(readxl);library(ggpubr);library(pheatmap)
##score
load("~/Downloads/scRNA/EC/EC/6_score/PCD_score_genesets.RData")
{
scRNA1 <- AddModuleScore(scRNA, features = markers[1], ctrl = 100, name = "Apoptosis")
scRNA1 <- AddModuleScore(scRNA1, features = markers[2], ctrl = 100, name = "Pyroptosis")
scRNA1 <- AddModuleScore(scRNA1, features = markers[3], ctrl = 100, name = "Ferroptosis")
scRNA1 <- AddModuleScore(scRNA1, features = markers[4], ctrl = 100, name = "Autophagy")
scRNA1 <- AddModuleScore(scRNA1, features = markers[5], ctrl = 100, name = "Necroptosis")
scRNA1 <- AddModuleScore(scRNA1, features = markers[6], ctrl = 100, name = "Parthanatos")
scRNA1 <- AddModuleScore(scRNA1, features = markers[7], ctrl = 100, name = "Entotic")
scRNA1 <- AddModuleScore(scRNA1, features = markers[8], ctrl = 100, name = "Netotic")
scRNA1 <- AddModuleScore(scRNA1, features = markers[9], ctrl = 100, name = "Cuproptosis")
scRNA1 <- AddModuleScore(scRNA1, features = markers[10], ctrl = 100, name = "Lysome_dependent")
scRNA1 <- AddModuleScore(scRNA1, features = markers[11], ctrl = 100, name = "Alkaliptosis")
scRNA1 <- AddModuleScore(scRNA1, features = markers[12], ctrl = 100, name = "Oxeiprosis")
}

##Organ
{
  score_ID <- data.frame(scoreID=paste(names(markers), "1", sep = ""))
  score <- FetchData(scRNA1, vars = c("Organ",score_ID$scoreID))
  score_mean <- aggregate(score[,2:13], by=list(type=score$Organ),mean)
  rownames(score_mean) <- score_mean$type
  score_mean <- score_mean[,-1]
  score_mean <- na.omit(score_mean)
  mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
  p1 = pheatmap(score_mean,scale = "none", border_color= "grey", number_color= "black",
                fontsize_number=14,fontsize_row=8,fontsize_col=9,cellwidth=14,
                cellheight=14,cluster_rows=T,cluster_cols=T,treeheight_row = 10,treeheight_col = 10,
                color= mycol,show_rownames=T)
  ggsave(filename = "PFC_score_Organ.pdf", p1,device = 'pdf', width = 10, height = 20, units = 'cm')
}

##
library(ggpubr)
table(scRNA1$Group,scRNA1$Organ)
table(scRNA1$Group)
scRNA1 = scRNA1[,scRNA1$Organ %in% c("Aorta","Bladder","Brain","Heart","Intestines","Kidney","Liver","Lung","Muscle","Pancreas","PIT","Retina","SAT","SC","Skin","Testis","Thyroid","Uterus","VAT")]
scRNA1$Organ <- as.factor(as.character(scRNA1$Organ))
table(scRNA1$Group,scRNA1$Organ)
scRNA1$Organ_Age <- paste(scRNA1$Organ, scRNA1$Age, sep = "_")
asd=data.frame(table(scRNA1$Organ_Age))
write.csv(asd,"Organ_Age.csv", row.names = F) 
cluster_celltype <-  data.frame(readxl::read_xlsx("Organ_Age.xlsx"))
current.cluster.ids <- cluster_celltype$Var1
new.cluster.ids <- cluster_celltype$Group
scRNA1$Group <- plyr::mapvalues(x = as.character(scRNA1$Organ_Age), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA1$Group,scRNA1$Organ)
scRNA1 = scRNA1[,scRNA1$Group %in% c("Adult","Aging")]
scRNA1$Group <- as.factor(as.character(scRNA1$Group))
table(scRNA1$Group,scRNA1$Organ)

##P.value
{
{
  #
  data <- FetchData(scRNA1, vars = c("Group","Apoptosis1","Organ"))
  p_ME <- compare_means(Apoptosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  Apoptosis <- data.frame(Organ=p_ME$Organ,Apoptosis=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Pyroptosis1","Organ"))
  p_ME <- compare_means(Pyroptosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  Pyroptosis <- data.frame(Organ=p_ME$Organ,Pyroptosis=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Ferroptosis1","Organ"))
  p_ME <- compare_means(Ferroptosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  Ferroptosis <- data.frame(Organ=p_ME$Organ,Ferroptosis=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Autophagy1","Organ"))
  p_ME <- compare_means(Autophagy1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  Autophagy <- data.frame(Organ=p_ME$Organ,Autophagy=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Necroptosis1","Organ"))
  p_ME <- compare_means(Necroptosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  Necroptosis <- data.frame(Organ=p_ME$Organ,Necroptosis=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Cuproptosis1","Organ"))
  p_ME <- compare_means(Cuproptosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  Cuproptosis <- data.frame(Organ=p_ME$Organ,Cuproptosis=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Parthanatos1","Organ"))
  p_ME <- compare_means(Parthanatos1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  Parthanatos <- data.frame(Organ=p_ME$Organ,Parthanatos=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Entotic1","Organ"))
  p_ME <- compare_means(Entotic1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  Entotic <- data.frame(Organ=p_ME$Organ,Entotic=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Netotic1","Organ"))
  p_ME <- compare_means(Netotic1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  Netotic <- data.frame(Organ=p_ME$Organ,Netotic=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Lysome_dependent1","Organ"))
  p_ME <- compare_means(Lysome_dependent1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  Lysome_dependent <- data.frame(Organ=p_ME$Organ,Lysome_dependent=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Alkaliptosis1","Organ"))
  p_ME <- compare_means(Alkaliptosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  Alkaliptosis <- data.frame(Organ=p_ME$Organ,Alkaliptosis1=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Oxeiprosis1","Organ"))
  p_ME <- compare_means(Oxeiprosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
  Oxeiprosis <- data.frame(Organ=p_ME$Organ,Oxeiprosis=p_ME$p.adj)
}
{
  p_ME_all <- merge(Apoptosis,  Pyroptosis, by="Organ",all=TRUE)
  p_ME_all <- merge(p_ME_all, Ferroptosis, by="Organ",all=TRUE)
  p_ME_all <- merge(p_ME_all, Autophagy, by="Organ",all=TRUE)
  p_ME_all <- merge(p_ME_all, Necroptosis, by="Organ",all=TRUE)
  p_ME_all <- merge(p_ME_all, Cuproptosis, by="Organ",all=TRUE)
  p_ME_all <- merge(p_ME_all, Parthanatos, by="Organ",all=TRUE)
  p_ME_all <- merge(p_ME_all, Entotic, by="Organ",all=TRUE)
  p_ME_all <- merge(p_ME_all, Netotic, by="Organ",all=TRUE)
  p_ME_all <- merge(p_ME_all, Lysome_dependent, by="Organ",all=TRUE)
  p_ME_all <- merge(p_ME_all, Alkaliptosis, by="Organ",all=TRUE)
  p_ME_all <- merge(p_ME_all, Oxeiprosis, by="Organ",all=TRUE)
}
head(p_ME_all)
write.csv(p_ME_all,"p_score_all.csv", na="0",row.names = F)
}
##-logP.value
{
  {
    #
    data <- FetchData(scRNA1, vars = c("Group","Apoptosis1","Organ"))
    p_ME <- compare_means(Apoptosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    Apoptosis <- data.frame(Organ=p_ME$Organ,Apoptosis=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Pyroptosis1","Organ"))
    p_ME <- compare_means(Pyroptosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    Pyroptosis <- data.frame(Organ=p_ME$Organ,Pyroptosis=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Ferroptosis1","Organ"))
    p_ME <- compare_means(Ferroptosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    Ferroptosis <- data.frame(Organ=p_ME$Organ,Ferroptosis=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Autophagy1","Organ"))
    p_ME <- compare_means(Autophagy1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    Autophagy <- data.frame(Organ=p_ME$Organ,Autophagy=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Necroptosis1","Organ"))
    p_ME <- compare_means(Necroptosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    Necroptosis <- data.frame(Organ=p_ME$Organ,Necroptosis=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Cuproptosis1","Organ"))
    p_ME <- compare_means(Cuproptosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    Cuproptosis <- data.frame(Organ=p_ME$Organ,Cuproptosis=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Parthanatos1","Organ"))
    p_ME <- compare_means(Parthanatos1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    Parthanatos <- data.frame(Organ=p_ME$Organ,Parthanatos=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Entotic1","Organ"))
    p_ME <- compare_means(Entotic1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    Entotic <- data.frame(Organ=p_ME$Organ,Entotic=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Netotic1","Organ"))
    p_ME <- compare_means(Netotic1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    Netotic <- data.frame(Organ=p_ME$Organ,Netotic=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Lysome_dependent1","Organ"))
    p_ME <- compare_means(Lysome_dependent1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    Lysome_dependent <- data.frame(Organ=p_ME$Organ,Lysome_dependent=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Alkaliptosis1","Organ"))
    p_ME <- compare_means(Alkaliptosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    Alkaliptosis <- data.frame(Organ=p_ME$Organ,Alkaliptosis1=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Oxeiprosis1","Organ"))
    p_ME <- compare_means(Oxeiprosis1~Group, data, group.by = "Organ",method = "wilcox.test", paired = FALSE)
    Oxeiprosis <- data.frame(Organ=p_ME$Organ,Oxeiprosis=-log10(p_ME$p.adj))
  }
  {
    p_ME_all <- merge(Apoptosis,  Pyroptosis, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, Ferroptosis, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, Autophagy, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, Necroptosis, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, Cuproptosis, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, Parthanatos, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, Entotic, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, Netotic, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, Lysome_dependent, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, Alkaliptosis, by="Organ",all=TRUE)
    p_ME_all <- merge(p_ME_all, Oxeiprosis, by="Organ",all=TRUE)
  }
  head(p_ME_all)
  write.csv(p_ME_all,"logp_score_all.csv", na="0",row.names = F)
}

##up_down
{
scRNA1$Organ_Group <- paste(scRNA1$Organ, scRNA1$Group, sep = "_")
score_ID <- data.frame(scoreID=paste(names(markers), "1", sep = ""))
score<-FetchData(scRNA1, vars = c("Organ_Group",score_ID$scoreID))
score_mean <- aggregate(score[,2:13], by=list(type=score$Organ_Group),mean)
rownames(score_mean) <- score_mean$type
score_mean <- score_mean[,-1]
score_mean <- na.omit(score_mean)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap(t(score_mean), show_colnames = T,scale = "row", cluster_cols = F,treeheight_col = 10, treeheight_row=15,color= mycol)
}


head(score_dif)
write.csv(score_dif,"score_all_diff.csv", na="0") #保存结果


library(pheatmap)
r2 <- read.csv('logp_score_all.csv',header=TRUE,row.names=1)
p2 <- read.csv('p_score_all.csv',header=TRUE,row.names=1)
p3 <- read.csv('scFEA_all_diff.csv',header=TRUE,row.names=1)

p2[p2 < 0.001] <- "***"
p2[p2 >= 0.001 & p2 < 0.01] <- "**"
p2[p2 >= 0.01 & p2 < 0.05] <- "*"
p2[p2 > 0.05] <- ""

#
p3[p3 > 0] <- 1
p3[p3 < 0] <- -1

#
r2 <- as.matrix(r2)
r2[is.infinite(r2)] <- 30
r2[r2 >= 30] <- 30
r2 <- r2*p3
#mycol<-colorRampPalette(c( "#0f86a9", "white", "#ed8b10"))(200)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)

#
p1 = pheatmap(r2,scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=12,fontsize_row=8,fontsize_col=9,cellwidth=15,
              cellheight=15,cluster_rows=T,cluster_cols=T,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "PCD1.pdf", p1,device = 'pdf', width = 15, height = 15, units = 'cm')
p1 = pheatmap(r2,scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=12,fontsize_row=8,fontsize_col=9,cellwidth=16,
              cellheight=16,cluster_rows=T,cluster_cols=F,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "PCD2.pdf", p1,device = 'pdf', width = 15, height = 15, units = 'cm')

################################################################################
# Extended data Fig. 10D
################################################################################
#scTenifoldKnk
load('Brain_Endo')
table(scRNA$Group)
scRNA_seq <- scRNA[,scRNA$Group %in% c("Adult")]
scRNA_seq$Group <- as.factor(as.character(scRNA_seq$Group))
X <- data.frame(scRNA_seq@assays$RNA@counts)
save(X ,file="Brain_Adult_KO_counts.RData")
##
library(scTenifoldKnk);library(Seurat);library(Matrix)
args = commandArgs(trailingOnly=TRUE)
scQC <- function(X, mtThreshold = 0.1, minLSize = 1000){
  require(Matrix)
  if(class(X) == 'Seurat'){
    countMatrix <- X@assays$RNA@counts
  } else {
    countMatrix <- X
  }
  librarySize <- colSums(countMatrix)
  countMatrix <- countMatrix[,librarySize >= minLSize]
  librarySize <- colSums(countMatrix)
  mtGenes <- grep('^MT-',toupper(rownames(countMatrix)))
  nGenes <- colSums(countMatrix != 0)
  
  genesLM <- lm(nGenes~librarySize)
  genesLM <- as.data.frame(predict(genesLM, data.frame(librarySize), interval = 'prediction'))
  
  if(isTRUE(length(mtGenes) > 0)){
    mtCounts <- colSums(countMatrix[grep('^MT-',toupper(rownames(countMatrix))),])
    mtProportion <- mtCounts/librarySize
    mtLM <- lm(mtCounts~librarySize)
    mtLM <- as.data.frame(predict(mtLM, data.frame(librarySize), interval = 'prediction'))
    selectedCells <- ((mtCounts > mtLM$lwr) & (mtCounts < mtLM$upr) & (nGenes > genesLM$lwr) & (nGenes < genesLM$upr) & (mtProportion <= mtThreshold) & (librarySize < 2 * mean(librarySize)))
  } else {
    selectedCells <- ((nGenes > genesLM$lwr) & (nGenes < genesLM$upr) & (librarySize < 2 * mean(librarySize)))
  }
  selectedCells <- colnames(countMatrix)[selectedCells]
  if(class(X) == 'Seurat'){
    X <- subset(X, cells = selectedCells)
  } else {
    X <- countMatrix[,selectedCells]
  }
  return(X)
}

# Running scTenifoldKnk
load("Brain_Adult_KO_counts.RData")
X <- scQC(X)
tempX <- X[,sample(seq_len(ncol(X)),2000)]
tempX <- tempX[rowMeans(X != 0) > 0.05,]
oX <- scTenifoldKnk(countMatrix = tempX, gKO = 'HIF3A',qc = TRUE, nCores = 16)
save(oX ,file="Brain_HIF3A_counts.RData")
write.csv(oX$diffRegulation, file = paste0('HIF3A_counts_Stability_',args[1], '.csv'), row.names = FALSE)
##
library(fgsea);library(ggplot2);library(enrichR);library(igraph)
source('/Users/niuruize/Downloads/scRNA/R_code/plotDR.R')
source('/Users/niuruize/Downloads/scRNA/R_code/plotKO.R')
source('/Users/niuruize/Downloads/scRNA/R_code/hsa2mmu_SYMBOL.R')

load('Brain_HIF3A_counts.RData')
Brain_HIF3A = oX
write.csv(Brain_HIF3A$diffRegulation, 'Brain_HIF3A_counts.csv')
dGenes <- Brain_HIF3A$diffRegulation$gene[Brain_HIF3A$diffRegulation$p.adj < 0.05]

#png('dr2_Brain_HIF3A.png', width = 2000, height = 2000, res = 300)
pdf("dr2_Brain_HIF3A.pdf", width = 5,height = 5)
plotDR(Brain_HIF3A)
dev.off()

#png('ego2_Brain_HIF3A.png', width = 3000, height = 3000, res = 300, bg = NA)
pdf("ego2_Brain_HIF3A.pdf", width = 10,height = 10)
X <- Brain_HIF3A
gKO <- 'HIF3A'
q <- 0.995
gList <- unique(c(gKO, X$diffRegulation$gene[X$diffRegulation$p.adj < 0.05]))
sCluster <- as.matrix(X$tensorNetworks$WT[gList,gList])
koInfo <- sCluster[gKO,]
sCluster[abs(sCluster) <= quantile(abs(sCluster), q)] <- 0
sCluster[gKO,] <- koInfo
diag(sCluster) <- 0
sCluster <-  reshape2::melt(as.matrix(sCluster))
colnames(sCluster) <- c('from', 'to', 'W')
sCluster <- sCluster[sCluster$W != 0,]
netPlot <- graph_from_data_frame(sCluster, directed = TRUE)
dPlot <- centr_degree(netPlot)$res
W <- rep(1,nrow(sCluster))
sG   <- (names(V(netPlot))[dPlot > 1])[-1]
W[sCluster$from %in% sG] <- 0.2
W[sCluster$to %in% sG] <- 0.2
W[sCluster$from %in% gKO] <- 1
W[sCluster$from %in% gKO & sCluster$to %in% sG] <- 0.8
set.seed(1)
layPlot <- layout_with_fr(netPlot, weights = W)
dPlot <- (dPlot/max(dPlot))*20
dbs <- listEnrichrDbs()
E <- enrichr(gList, c("GO_Biological_Process_2023","KEGG_2021_Human"))
E <- do.call(rbind.data.frame, E)
#E <- E[E$Adjusted.P.value < 0.05,]
E <- E[order(E$Adjusted.P.value),]
E$Term <- unlist(lapply(strsplit(E$Term,''), function(X){
  X[1] <- toupper(X[1])
  X <- paste0(X,collapse = '')
  X <- gsub('\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+','',X)
  X <- gsub("'s",'',X)
  X <- unlist(strsplit(X,','))[1]
  X <- gsub('[[:blank:]]$','',X)
  return(X)
}))
write.csv(E,"HIF3A_EGO.csv") #保存结果
E <- E[E$Term %in% c('Sequestering Of Actin Monomers','Cardiac Cell Development','Cellular Response To Metal Ion','','Regulation Of Blood Vessel Endothelial Cell Migration','Vasoconstriction','Regulation Of Endothelial Cell Chemotaxis'),]
#E <- E[c(1,2,3,4),]
tPlot <- strsplit(E$Genes, ';')
pPlot <- matrix(0,nrow = length(V(netPlot)), ncol = nrow(E))
rownames(pPlot) <- toupper(names(V(netPlot)))
for(i in seq_along(tPlot)){
  pPlot[unlist(tPlot[i]),i] <- 1
}
pPlot <- lapply(seq_len(nrow(pPlot)), function(X){as.vector(pPlot[X,])})
names(pPlot) <- names(V(netPlot))
tPlot <- unique(unlist(tPlot))
eGenes <- toupper(names(V(netPlot))) %in% tPlot
vColor <- rgb(195/255, 199/255, 198/255 ,0.3)
pieColors <- list(hcl.colors(nrow(E), palette = 'Zissou 1', alpha = 0.7))
par(mar=c(4,0,0,0), xpd = TRUE)
suppressWarnings(plot(netPlot,
                      layout = layPlot,
                      edge.arrow.size=.2,
                      vertex.label.color="black",
                      vertex.shape = ifelse(eGenes,'pie','circle'),
                      vertex.pie = pPlot,
                      vertex.size = 10+dPlot,
                      vertex.pie.color=pieColors,
                      vertex.label.family="Times",
                      vertex.label.font=ifelse(eGenes,2,1),
                      edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue'),
                      edge.curved = ifelse(W == 0.2, 0, 0.1),
                      vertex.color = vColor,
                      vertex.frame.color = NA))
sigLevel <- formatC(E$Adjusted.P.value, digits = 2, format = 'g', width = 0, drop0trailing = TRUE)
gSetNames <- lengths(strsplit(E$Genes, ';'))
gSetNames <- paste0('(', gSetNames,') ', E$Term, ' FDR = ', sigLevel)
legend(x = -1.05, y = -1.05, legend = gSetNames, bty = 'n', ncol = 2, cex = 1, col = unlist(pieColors), pch = 16)
dev.off()


zBrain_HIF3A <- Brain_HIF3A$diffRegulation$Z
names(zBrain_HIF3A) <- toupper(Brain_HIF3A$diffRegulation$gene)

MGI <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human')
set.seed(1)
E <- fgseaMultilevel(MGI, zBrain_HIF3A)

#png('gsea1_Brain_HIF3A.png', width = 1000, height = 1000, res = 300)
pdf("gsea1_Brain_HIF3A.pdf", width = 5,height = 5)
gSet <- 'Fluid shear stress and atherosclerosis'
plotEnrichment(MGI[[gSet]], zBrain_HIF3A) +
  labs(
    title = 'Fluid shear stress and atherosclerosis',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea2_Brain_HIF3A.png', width = 1000, height = 1000, res = 300)
pdf("gsea2_Brain_HIF3A.pdf", width = 5,height = 5)
gSet <- 'Complement and coagulation cascades'
plotEnrichment(MGI[[gSet]], zBrain_HIF3A) +
  labs(
    title = 'Complement and coagulation cascades',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

write.csv(data.frame(lapply(E, as.character)),"HIF3A_GSEA.csv") 































