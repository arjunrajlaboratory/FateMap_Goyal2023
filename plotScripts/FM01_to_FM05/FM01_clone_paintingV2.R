########################################################################################
#Description: Clone painting on big cluster
#Author: Maalavika Pillai
#Version: 2
#Edited on: 11/16/22
#Includes distribution of the clone size per sample and how many were singletons and
#Scoring for cell cycle/ apoptosis and
#Expression of mutant genes
########################################################################################

library(ggplot2)
library(dplyr)
library(ggpubr)
library(svglite)
library(AUCell)

#Load data and barcodes
setwd("~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM01/Data/")
fm01 <- readRDS("s1s2ScTransform_50pcs_filter.rds")
barcodes <- read.table('barcodeCellID.tsv', stringsAsFactors=F, header = T, sep="\t")
barcodes$cellID <- paste0(barcodes$sampleNum, "_", barcodes$cellID)

#Get sizes of barcodes
barcodes <- barcodes %>%
  group_by(BC50StarcodeD8) %>%
  mutate(nSize = length(cellID))

#Distribution of clone sizes 
plot <- ggplot(barcodes) + 
  geom_histogram(aes(x=nSize), binwidth = 10,color="black", fill="black",boundary = 0)+
  theme_classic()+
  theme(text = element_text(size = 22))+ 
  ylab("Frequency")+
  xlab('Clone Size')+
  scale_x_continuous(breaks = seq(0,160,10))+
  annotate("text", x=80, y=2000, label= paste0("Total number of barcoded samples = ", nrow(barcodes)), size = 6) + 
  annotate("text", x=80, y=1930, label= paste0("Total number of singletons = ", sum(barcodes$nSize==1), "  (",round(100*sum(barcodes$nSize==1)/ nrow(barcodes),2),"%)"), size=6)+
  annotate("text", x=80, y=1860, label= paste0("Largest clone fraction = ", round(max(barcodes$nSize)/ nrow(barcodes),4) ), size=6)
ggsave(plot = plot, "../Plots/CloneSizeHistogram.svg", width = 10, height = 10)

fm01_barcoded <- fm01[,colnames(fm01) %in% barcodes$cellID]
rownames(barcodes) <- barcodes$cellID
barcodes <- barcodes[colnames(fm01_barcoded),]
fm01_barcoded$barcodes <- barcodes$BC50StarcodeD8
fm01_barcoded$nSize <- barcodes$nSize

#Plot mutant genes
fm01 <- FindClusters(fm01, resolution = 0.04)
plot <- DoHeatmap(fm01, features = c("KRT18", "GXYLT1", "ANKRD36", "NBPF10", "PDE4DIP"))
ggsave(plot = plot, "../Plots/MutantGenes_2clusters.svg", width = 10, height = 10)

plot <-DoHeatmap(fm01, features = c("KRT18", "GXYLT1", "ANKRD36", "NBPF10", "PDE4DIP"), group.by = "integrated_snn_res.0.6")
ggsave(plot = plot, "../Plots/MutantGenes_15clusters.svg", width = 10, height = 10)

plot <- FeaturePlot(fm01, features = c("KRT18", "GXYLT1", "ANKRD36", "NBPF10", "PDE4DIP", "LEF1", "FMN1","HLA-A", "SERPINE1"), ncol = 2)
ggsave(plot = plot, "../Plots/MutantGenes_UMAP_AllCells.svg", width = 10, height = 25)

plot <- FeaturePlot(fm01_barcoded, features = c("KRT18", "GXYLT1", "ANKRD36", "NBPF10", "PDE4DIP"))
ggsave(plot = plot, "../Plots/MutantGenes_UMAP_barcodedCells.svg", width = 10, height = 15)


#Distribution of clone sizes (for)
barcodes_S2 <-read.table('barcodeCellID_S2.tsv', stringsAsFactors=F, header = T, sep="\t")
barcodes_S2$cellID <- paste0(barcodes_S2$sampleNum, "_", barcodes_S2$cellID)

#Get sizes of barcodes
barcodes_S2 <- barcodes_S2 %>%
  group_by(BC50StarcodeD8) %>%
  mutate(nSize = length(cellID))

plot <- ggplot(barcodes_S2) + 
  geom_histogram(aes(x=nSize), binwidth = 10,color="black", fill="black",boundary = 0)+
  theme_classic()+
  theme(text = element_text(size = 22))+ 
  ylab("Frequency")+
  xlab('Clone Size')+
  scale_x_continuous(breaks = seq(0,160,10))+
  annotate("text", x=80, y=2000, label= paste0("Total number of barcoded samples = ", nrow(barcodes_S2)), size = 6) + 
  annotate("text", x=80, y=1930, label= paste0("Total number of singletons = ", sum(barcodes_S2$nSize==1), "  (",round(100*sum(barcodes_S2$nSize==1)/ nrow(barcodes_S2),2),"%)"), size=6)+
  annotate("text", x=80, y=1860, label= paste0("Largest clone fraction = ", round(max(barcodes_S2$nSize)/ nrow(barcodes_S2),4) ), size=6)
ggsave(plot = plot, "../Plots/CloneSizeHistogram_S2.svg", width = 10, height = 10)

#Calculate cell cycle scores and apoptosis scores using AUCell
DefaultAssay(fm01_barcoded) <- "RNA"
exprMatrix <- as.matrix(GetAssayData(object =fm01_barcoded, slot = "counts")) #Extract counts
cells_rankings <- AUCell_buildRankings(exprMatrix)

#Read genesets
geneSets <- list()
genelist <- list.files("CellCycle_Apoptosis_Genes/")
for (i in genelist) {
  geneSets[[which(genelist==i)]] <- as.character(read.delim(paste0("CellCycle_Apoptosis_Genes/",i),sep = "\t", header = F)[,-c(1:2)])
}
names(geneSets) <- gsub(".v.*", "", genelist)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
AUCscores <- getAUC(cells_AUC) #Extract scores for each sample
fm01_barcoded@meta.data <- cbind(fm01_barcoded@meta.data, as.data.frame(t(AUCscores))) #Add to metadat of Seurat object

DefaultAssay(fm01_barcoded) <- "integrated"
#UMAP plot for various scores
plot <- FeaturePlot(fm01_barcoded, gsub(".v.*", "", genelist))
ggsave(plot = plot, "../Plots/CellCycleApoptosis_AUCell.svg", width = 10, height = 10)

#UMAP for singlets
fm01_barcoded$Singlet <- fm01_barcoded$nSize ==1
DimPlot(fm01_barcoded, group.by = 'Singlet')
DimPlot(fm01_barcoded)

#Select larger cluster
fm01_barcoded <- FindNeighbors(fm01_barcoded)
fm01_barcoded <- FindClusters(fm01_barcoded, resolution = 0.04) #2 clusters
DimPlot(fm01_barcoded)
fm01_barcoded_large <- fm01_barcoded[,fm01_barcoded$seurat_clusters==0&!fm01_barcoded$Singlet]
fm01_barcoded_large <- RunPCA(fm01_barcoded_large)
fm01_barcoded_large <- RunUMAP(fm01_barcoded_large, dims = 1:50)
fm01_barcoded_large <- FindNeighbors(fm01_barcoded_large)
fm01_barcoded_large <- FindClusters(fm01_barcoded_large)
plot <- DimPlot(fm01_barcoded_large)
ggsave(plot = plot, "../Plots/LargeCluster_subcluster.svg", width = 10, height = 10)

fm01_barcoded_large$barcodes <- factor(fm01_barcoded_large$barcodes, levels = unique(fm01_barcoded_large$barcodes), labels = 1:length(unique(fm01_barcoded_large$barcodes)))
plot <- DimPlot(fm01_barcoded_large, group.by =  "barcodes") +NoLegend()
ggsave(plot = plot, "../Plots/LargeCluster_AllBarcodes.svg", width = 10, height = 10)

#Select 3 random clusters
set.seed(100)
samp <- sample(levels(fm01_barcoded_large$barcodes), 25)
plots <- list()
for (i in samp) {
  cell_names = colnames(fm01_barcoded_large)[fm01_barcoded_large$barcodes==i]
  plots[[which(samp ==i)]] <- DimPlot(fm01_barcoded_large, cells.highlight = cell_names, cols.highlight = "hotpink3") +NoLegend()
}
plot <- ggarrange(plotlist = plots, nrow=5, ncol =5)
ggsave(plot = plot, "../Plots/LargeCluster_25RandomBarcodes.svg", width = 15, height = 15)
