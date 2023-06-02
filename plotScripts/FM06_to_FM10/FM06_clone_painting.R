########################################################################################
#Description: Clone painting of FM06
#Author:  Maalavika Pillai
#Version: 1
#Edited on: 10/28/22
########################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(svglite)
#Load data and barcodes
setwd("~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2021_FM06/extracted/FM06")
fm06 <- readRDS("s1s2Naive_ScTransform_50pcs_filter.rds")
barcodes <- read.table('../FM06_barcode/barcodeCellID.tsv', stringsAsFactors=F, header = T, sep="\t")
barcodes$cellID <- paste0(barcodes$sampleNum, "_", barcodes$cellID)

#Get sizes of barcodes
barcodes <- barcodes %>%
  group_by(BC50StarcodeD8) %>%
  mutate(nSize = length(cellID))

fm06_barcoded <- fm06[,colnames(fm06) %in% barcodes$cellID]
rownames(barcodes) <- barcodes$cellID
barcodes <- barcodes[colnames(fm06_barcoded),]
fm06_barcoded$barcodes <- barcodes$BC50StarcodeD8
fm06_barcoded$nSize <- barcodes$nSize

#UMAP for singlets
fm06_barcoded$Singlet <- fm06_barcoded$nSize ==1
fm06_barcoded$Singlet <- factor(fm06_barcoded$Singlet , levels = c(TRUE, FALSE)) 
plot <- DimPlot(fm06_barcoded, group.by = 'Singlet' )+ 
  scale_color_manual(values= c( "hotpink3", "lightgrey")) + 
  NoLegend()
ggsave(plot = plot, filename = paste0("../../plot/FM06_UMAP_singlet.svg"), width = 10, height = 10)

fm06_barcoded_notSingleton <-  fm06_barcoded[,fm06_barcoded$Singlet==FALSE ]
Idents(object = fm06_barcoded_notSingleton) <-fm06_barcoded_notSingleton$barcodes
cloneMarkers <- FindAllMarkers(fm06_barcoded_notSingleton, only.pos = T, min.pct = 0.5, logfc.threshold = 1)
cloneMarkers <- cloneMarkers[cloneMarkers$avg_log2FC>1 & cloneMarkers$p_val < 0.05,]
write.table(cloneMarkers, paste0("CloneGenes_FC1_p0.05_AllSinglets.tsv"), sep = "\t", quote = F)

#UMAP for clusters
plot <- DimPlot(fm06_barcoded)
ggsave(plot = plot, filename = paste0("../../plot/FM06_UMAP_seuratCluster.svg"), width = 10, height = 10)

#Gene painting
genes <- c("ACTA2","IFIT2", "VCAM1", "NGFR", "MITF", "AXL")
DefaultAssay(fm06_barcoded) <- "SCT"
plot <- FeaturePlot(fm06_barcoded, features = genes, slot = "data" )
ggsave(plot = plot, filename = paste0("../../plot/FM06_FeaturePlot_genePainting.svg"), width = 7, height = 10)

#UMAP for clones present in both samples/twins
barcode_twins <- vector()

for (i in unique(barcodes$BC50StarcodeD8)) {
  if(all(c("S1","S2") %in% barcodes$sampleNum[barcodes$BC50StarcodeD8 ==i ] )){
    barcode_twins <- c(barcode_twins, i) 
  }
}

fm06_barcoded$Sample <- substr(colnames(fm06_barcoded),1,2)
set.seed(100)
randSampBarcodes <- sample(barcode_twins, 25) 
plots <- list()
for (i in randSampBarcodes) {
  cell_names_S1  <-colnames(fm06_barcoded[,fm06_barcoded$barcodes == i & fm06_barcoded$Sample=="S1"])
  cell_names_S2  <-colnames(fm06_barcoded[,fm06_barcoded$barcodes == i & fm06_barcoded$Sample=="S2"])
  plots[[which(randSampBarcodes == i)]] <- DimPlot(fm06_barcoded, cells.highlight = list(cell_names_S1, cell_names_S2), 
                                                  cols.highlight = c("coral2", "darkmagenta")) + ggtitle(i) +NoLegend()
} 
plot <- ggarrange(plotlist = plots, nrow = 5, ncol =5)
ggsave(plot = plot, filename = paste0("../../plot/FM06_Twin_memory.svg"), width = 20, height = 20)

#Specific examples
set.seed(100)
randSampBarcodes <- sample(barcode_twins, 25) [c(24,19,14,1,2)]
plots <- list()
for (i in randSampBarcodes) {
  cell_names_S1  <-colnames(fm06_barcoded[,fm06_barcoded$barcodes == i & fm06_barcoded$Sample=="S1"])
  cell_names_S2  <-colnames(fm06_barcoded[,fm06_barcoded$barcodes == i & fm06_barcoded$Sample=="S2"])
  plots[[which(randSampBarcodes == i)]] <- DimPlot(fm06_barcoded, cells.highlight = list(cell_names_S1, cell_names_S2), 
                                                   cols.highlight = c("coral2", "darkmagenta")) + ggtitle(i) +NoLegend()
} 
plot <- ggarrange(plotlist = plots, nrow = 1)
ggsave(plot = plot, filename = paste0("../../plot/FM06_Twin_memory_specificClones.svg"), width = 30, height =8)



#select barcodes with >10 cells
fm06_rand <- fm06_barcoded[,fm06_barcoded$nSize > 10]
plots <- list()
for (bc_samp in unique(fm06_rand$barcodes)) {
  cell_names <- colnames(fm06_rand[,fm06_rand$barcodes == bc_samp])
  plots[[which(unique(fm06_rand$barcodes) == bc_samp)]] <- DimPlot(fm06, cells.highlight = cell_names  )+ NoLegend()
}
plot <- ggarrange(plotlist = plots,  ncol=5, nrow = 5)
ggsave(plot = plot, filename = paste0("../../plot/FM06_UMAP_cloneSize_10+.svg"), width = 10, height = 10)


#Barcodes for 3-10 clone siz

barcodes <- list()
for(i in 3:10){
  fm06_rand <- fm06_barcoded[,fm06_barcoded$nSize == i]
  plots <- list()
  barcodes_size <- unique(fm06_rand$barcodes)
  n <- length(barcodes_size)
  n <- (floor(n/5))
  n <- min(n,5)
  set.seed(100)
  barcodes_size <- sample(barcodes_size, (n*5))
  barcodes[[i-2]] <- barcodes_size
  for (bc_samp in barcodes_size) {
    cell_names <- colnames(fm06_rand[,fm06_rand$barcodes == bc_samp])
    plots[[which(barcodes_size == bc_samp)]] <- DimPlot(fm06, cells.highlight = cell_names  ) + NoLegend()
  }
  plot <- ggarrange(plotlist = plots,  ncol=5, nrow = n) + ggtitle(paste0("Clone size = ", i))
  ggsave(plot = plot, filename = paste0("../../plot/FM06_UMAP_cloneSize", i, ".svg"), width = 10, height = 10)
}


#Highlighting 5 random barcodes
barcodes[[9]] <-  unique(fm06_barcoded$barcodes[fm06_barcoded$nSize > 10])

plots <- list()
barcode <- barcodes[[2]][12]
cell_names <- colnames(fm06_barcoded[,fm06_barcoded$barcodes == barcode])
plots[[1]] <- DimPlot(fm06, cells.highlight = cell_names , cols.highlight = "hotpink3") + NoLegend() + ggtitle(barcode)


barcode <- barcodes[[9]][13]
cell_names <- colnames(fm06_barcoded[,fm06_barcoded$barcodes == barcode])
plots[[2]] <- DimPlot(fm06, cells.highlight = cell_names , cols.highlight = "hotpink3" ) + NoLegend()+ ggtitle(barcode)

barcode <- barcodes[[4]][16]
cell_names <- colnames(fm06_barcoded[,fm06_barcoded$barcodes == barcode])
plots[[3]] <- DimPlot(fm06, cells.highlight = cell_names , cols.highlight = "hotpink3" ) + NoLegend()+ ggtitle(barcode)

barcode <- barcodes[[5]][9]
cell_names <- colnames(fm06_barcoded[,fm06_barcoded$barcodes == barcode])
plots[[4]] <- DimPlot(fm06, cells.highlight = cell_names , cols.highlight = "hotpink3" ) + NoLegend()+ ggtitle(barcode)

barcode <- barcodes[[6]][1]
cell_names <- colnames(fm06_barcoded[,fm06_barcoded$barcodes == barcode])
plots[[5]] <- DimPlot(fm06, cells.highlight = cell_names , cols.highlight = "hotpink3" ) + NoLegend()+ ggtitle(barcode)

plot <- ggarrange(plotlist = plots, nrow = 1)
ggsave(plot = plot, filename = paste0("../../plot/FM06_UMAP_specificClones.svg"), width = 30, height = 8)
