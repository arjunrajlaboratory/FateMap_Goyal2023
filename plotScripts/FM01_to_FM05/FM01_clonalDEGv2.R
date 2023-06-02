########################################################################################
#Description: Identify DEGs in clones
#Aim:Create  histogram for clone size, identify DEGs in top n clones, compare to DEGs
#from clusters
#Author: Maalavika Pillai
#Version: 2.1
#Edits: Assumes all singlets as a single barcode/clone
#Created on: 10/27/22
########################################################################################
library(ggplot2)
library(tidyverse)
library(reshape2)
library(svglite)
library(Seurat)
library(dplyr)
library(VennDiagram)
library(RColorBrewer)
library(svglite)

# Load the relevant datasets (filteredbarcode matrix files)

home1Directory <-"~/Library/CloudStorage/OneDrive-NorthwesternUniversity/Fatemap Revisions/Data/Drug treated and naive/FM01/"
dataDirectory <-"~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM01_DEGandVarGenes/Data/"
plotDirectory <-"~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM01_DEGandVarGenes/Plots/"

##FM01_sct_integratedData_sample1sample2
fm01_sample1_2 = readRDS(file = paste0(home1Directory,'s1s2ScTransform_50pcs_filter.rds'))  

#Barcodes
barcode50 = read.delim(file = paste0(home1Directory, 'barcodeCellID.tsv'))
barcode50$cellID <- paste0(barcode50$sampleNum,'_', barcode50$cellID)
barcode50 <- barcode50[barcode50$cellID %in% (colnames(fm01_sample1_2)),]


#Histogram of clone size
barcode50 <- barcode50 %>%
  group_by(BC50StarcodeD8) %>%
  mutate(nSize = length(cellID))
barcode50$Singlet <- barcode50$nSize ==1


#Order barcodes and select cells with clone size > 1
fm01_barcoded <- fm01_sample1_2[, colnames(fm01_sample1_2) %in% barcode50$cellID]
rownames(barcode50) <- barcode50$cellID
barcode50 <- barcode50[colnames(fm01_barcoded),]
fm01_barcoded$barcodes <- barcode50$BC50StarcodeD8
fm01_barcoded$barcodes[barcode50$Singlet] <- "Singlet" 

#DEG analysis
Idents(object = fm01_barcoded) <- fm01_barcoded$barcodes
cloneMarkers <- FindAllMarkers(fm01_barcoded, only.pos = T, min.pct = 0.5, logfc.threshold = 1)
cloneMarkers <- cloneMarkers[cloneMarkers$avg_log2FC>1 & cloneMarkers$p_val < 0.05,]
write.table(cloneMarkers, paste0(dataDirectory, "CloneGenes_FC1_p0.05_AllSinglets.tsv"), sep = "\t", quote = F)
cloneMarkers <- read.delim(paste0(dataDirectory, "CloneGenes_FC1_p0.05_AllSinglets.tsv"), sep = "\t")
#sum(cloneMarkers$cluster=="Singlet")
cloneMarkersUnique <- unique(cloneMarkers$gene)

#Read cluster DEGs 
clusterMarkers <- read.delim(paste0(dataDirectory, "scTransformMarkers_snn06.tsv"), sep="\t")
clusterMarkers <- clusterMarkers$gene[clusterMarkers$avg_logFC>1&clusterMarkers$p_val<0.05]
clusterMarkers <- unique(clusterMarkers)

#Read top 1000 variable genes
varGenes <- read.delim(paste0(dataDirectory, "topVarGenes_1000.tsv"), sep="\t")$x

#VennDiagram for all 3
myCol <- brewer.pal(3, "Pastel2")
# Chart
venn.diagram(
  x = list(varGenes, clusterMarkers, cloneMarkersUnique),
  category.names = c( "Variable genes ", "Cluster DEGs", "Clone DEGs" ),
  filename = paste0(plotDirectory, "VarGenes_Cluster_CloneDEGSinglet_Venn.png" ),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 550 , 
  width = 550 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-23, 23, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans"
)


venn.diagram(
  x = list(clusterMarkers, cloneMarkersUnique),
  category.names = c( "Cluster DEGs", "Clone DEGs" ),
  filename = paste0(plotDirectory, "Cluster_CloneDEGSinglet_Venn.png" ),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 550 , 
  width = 550 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-23, 23),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)


#Check change in clonal genes
cloneGenes1 <- read.delim(paste0(dataDirectory, "CloneGenes_FC1_p0.05.tsv"), sep = "\t")
cloneGenes2<- read.delim(paste0(dataDirectory, "CloneGenes_FC1_p0.05_AllSinglets.tsv"), sep = "\t")
venn.diagram(
  x = list(cloneGenes1$gene, cloneGenes2$gene),
  category.names = c( "Clone DEG Original", "Clone DEGs Singlet" ),
  filename = paste0(plotDirectory, "CloneGenes_Comparison_Singlet_Venn.png" ),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 550 , 
  width = 550 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-23, 23),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)
