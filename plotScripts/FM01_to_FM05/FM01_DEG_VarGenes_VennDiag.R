########################################################################################
#Description: Gene overlap between DEG and variable genes for FM01
#Aim: Venn diagram for variable genes and DEG
#Author: Maalavika Pillai
#Version: 1 , adaplted from  Correlation_clone_vs_random_all.R
#Created on: 10/25/22
########################################################################################

library(Seurat)
library(ggpubr)
library(dplyr)
library(rstatix)
library(ggplot2)
library(tidyr)
library(VennDiagram)
library(RColorBrewer)
library(svglite)

#Load data and barcodes
setwd("~/OneDrive - Northwestern University/Fatemap Revisions/Data/Drug treated and naive/FM01/")
fm01 <- readRDS("s1s2ScTransform_50pcs_filter.rds")
barcodes <- read.table('../barcodeCellID.tsv', stringsAsFactors=F, header = T, sep="\t")
barcodes$cellID <- paste0(barcodes$sampleNum, "_", barcodes$cellID)

#Select samples that have barcode information available
df <- as.data.frame(t(as.matrix(fm01@assays$integrated@scale.data)))

#Remove barcodes with single cellIDs and only those absent in df
barcodes <- barcodes %>%
  filter(cellID %in% rownames(df)) %>%
  group_by(BC50StarcodeD8) %>%
  filter(length(cellID) >= 2)
df <- df[barcodes$cellID,]
barcodes$BC50StarcodeD8 <- factor(barcodes$BC50StarcodeD8, labels = 1:length(unique(barcodes$BC50StarcodeD8)))

#For DEG analysis specifically 
DEG_main <-read.delim("scTransformMarkers_snn06.tsv", sep="\t")
cut_off <- c(2,0.01)
DEG <- DEG_main$gene[DEG_main$avg_logFC>cut_off[1] & DEG_main$p_val_adj <cut_off[2]]
write.table(DEG_main[DEG_main$gene %in% DEG,], paste0("DEG_FC", cut_off[1],"_pval",cut_off[2],".tsv"),sep = "\t", quote = F, row.names = F)


#For variable genes
nVarGenes = 1000
fm01 <- FindVariableFeatures(fm01, assay="RNA", nfeatures = nVarGenes)
topVarGenes <- VariableFeatures(fm01, assay="RNA")
write.table(topVarGenes, paste0("topVarGenes_", nVarGenes, ".tsv"),sep = "\t", quote = F, row.names = F)


#Venn Diagram
myCol <- brewer.pal(2, "Pastel2")
# Chart
venn.diagram(
  x = list(topVarGenes, DEG),
  category.names = c( "Variable genes ", "DEG" ),
  filename = paste0("DEG_FC", cut_off[1], "_pval",cut_off[2],"DEG_",nVarGenes,".png" ),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c('#CBD5E8', '#FDCDAC'),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-15, 15),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)
