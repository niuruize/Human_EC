################################################################################
# Extended data Fig. 4A
################################################################################
# The code for the visualization is shown in pyscenic.ipynb

################################################################################
# Supplementary Figure 4D
################################################################################
#scTenifoldKnk
load("Placenta_EC.RData")
scRNA_seq <- scRNA
X <- data.frame(scRNA_seq@assays$RNA@counts)
save(X ,file="Placenta_KO_counts.RData")
# KO
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
load("Placenta_KO_counts.RData")
X <- scQC(X)
tempX <- X[,sample(seq_len(ncol(X)),140)]
tempX <- tempX[rowMeans(X != 0) > 0.05,]
oX <- scTenifoldKnk(countMatrix = tempX, gKO = 'ARID3A',qc = FALSE, nCores = 16)
save(oX ,file="Placenta_ARID3A_data.RData")
write.csv(oX$diffRegulation, file = paste0('ARID3A_data_Stability_',args[1], '.csv'), row.names = FALSE)
##
library(fgsea);library(ggplot2);library(enrichR);library(igraph)
source('/Users/niuruize/Downloads/scRNA/R_code/plotDR.R')
source('/Users/niuruize/Downloads/scRNA/R_code/plotKO.R')
source('/Users/niuruize/Downloads/scRNA/R_code/hsa2mmu_SYMBOL.R')

load('Placenta_ARID3A_data.RData')
Placenta_ARID3A = oX
write.csv(Placenta_ARID3A$diffRegulation, 'Placenta_ARID3A_counts.csv')
dGenes <- Placenta_ARID3A$diffRegulation$gene[Placenta_ARID3A$diffRegulation$p.adj < 0.05]

#png('dr2_Placenta_ARID3A.png', width = 2000, height = 2000, res = 300)
pdf("dr2_Placenta_ARID3A.pdf", width = 5,height = 5)
plotDR(Placenta_ARID3A)
dev.off()

#png('ego2_Placenta_ARID3A.png', width = 3000, height = 3000, res = 300, bg = NA)
pdf("ego2_Placenta_ARID3A.pdf", width = 10,height = 10)
X <- Placenta_ARID3A
gKO <- 'ARID3A'
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
write.csv(E,"ARID3A_EGO.csv") 
E <- E[E$Term %in% c('Aorta Morphogenesis','Aorta Development','Blood Vessel Morphogenesis','Venous Blood Vessel Development',
                     'Vascular Wound Healing','Maintenance Of Blood-Brain Barrier','Sprouting Angiogenesis',
                     'Allograft rejection','Antigen processing and presentation','Hematopoietic cell lineage',
                     'Regulation Of Immune Effector Process','MHC Class II Protein Complex Assembly'),]
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

zPlacenta_ARID3A <- Placenta_ARID3A$diffRegulation$Z
names(zPlacenta_ARID3A) <- toupper(Placenta_ARID3A$diffRegulation$gene)

MGI <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human')
set.seed(1)
E <- fgseaMultilevel(MGI, zPlacenta_ARID3A)

#png('gsea1_Placenta_ARID3A.png', width = 1000, height = 1000, res = 300)
pdf("gsea1_Placenta_ARID3A.pdf", width = 5,height = 5)
gSet <- 'Graft-versus-host disease'
plotEnrichment(MGI[[gSet]], zPlacenta_ARID3A) +
  labs(
    title = 'Graft-versus-host disease',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea2_Placenta_ARID3A.png', width = 1000, height = 1000, res = 300)
pdf("gsea2_Placenta_ARID3A.pdf", width = 5,height = 5)
gSet <- 'Hematopoietic cell lineage'
plotEnrichment(MGI[[gSet]], zPlacenta_ARID3A) +
  labs(
    title = 'Hematopoietic cell lineage',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea2_Placenta_ARID3A.png', width = 1000, height = 1000, res = 300)
pdf("gsea3_Placenta_ARID3A.pdf", width = 5,height = 5)
gSet <- 'Cell adhesion molecules (CAMs)'
plotEnrichment(MGI[[gSet]], zPlacenta_ARID3A) +
  labs(
    title = 'Cell adhesion molecules (CAMs)',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea2_Placenta_ARID3A.png', width = 1000, height = 1000, res = 300)
pdf("gsea4_Placenta_ARID3A.pdf", width = 5,height = 5)
gSet <- 'Antigen processing and presentation'
plotEnrichment(MGI[[gSet]], zPlacenta_ARID3A) +
  labs(
    title = 'Antigen processing and presentation',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

write.csv(data.frame(lapply(E, as.character)),"ARID3A_GSEA.csv") 















