########################################################################################
#Description: Clone painting of FM10
#Author: Maalavika Pillai
#Version: 1
#Edited on: 10/28/22
########################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(svglite) 
library(ggrepel)

dataDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM10/jointAnalysis/seurat/"
plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM10/jointAnalysis/seuratOutput/"

#Load data and barcodes
fm10 <- readRDS(paste0(dataDirectory, "s1s2ScTransform_50pcs_filter.rds"))
barcodes <- read.table(paste0(dataDirectory,'barcodeCellID.tsv'), stringsAsFactors=F, header = T, sep="\t")
barcodes$cellID <- paste0(barcodes$sampleNum, "_", barcodes$cellID)

#Get sizes of barcodes
barcodes <- barcodes %>%
  group_by(BC50StarcodeD8) %>%
  mutate(nSize = length(cellID))

fm10 <- RenameCells(fm10, new.names = substr(colnames(fm10),1,nchar(colnames(fm10))-2))
fm10_barcoded <- fm10[,colnames(fm10) %in% barcodes$cellID]
rownames(barcodes) <- barcodes$cellID
barcodes <- barcodes[colnames(fm10_barcoded),]
fm10_barcoded$barcodes <- barcodes$BC50StarcodeD8
fm10_barcoded$nSize <- barcodes$nSize

#UMAP for singlets
fm10_barcoded$Singlet <- fm10_barcoded$nSize ==1
fm10_barcoded$Singlet <- factor(fm10_barcoded$Singlet , levels = c(TRUE, FALSE)) 
plot <- DimPlot(fm10_barcoded, group.by = 'Singlet' )+ 
  scale_color_manual(values= c( "hotpink3", "lightgrey")) + 
  NoLegend()
ggsave(plot = plot, filename = paste0(plotDirectory ,"FM10_UMAP_singlet.svg"), width = 10, height = 10)

#UMAP for clusters
plot <- DimPlot(fm10_barcoded)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM10_UMAP_samples.svg"), width = 10, height = 10)
fm10_barcoded <- FindNeighbors(fm10_barcoded)
fm10_barcoded <- FindClusters(object=fm10_barcoded, resolution = 0.6, verbose = FALSE)
plot <- DimPlot(fm10_barcoded)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM10_UMAP_seuratCluster.svg"), width = 10, height = 10)

#Gene markers for clusters
fm10_barcoded.markers <- FindAllMarkers(fm10_barcoded, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fm10_barcoded.markersSubset <- fm10_barcoded.markers %>% filter(p_val_adj <0.05)
list = c("ACTG2", "ACTA2", "MYOCD", "TAGLN", "IFIT2", "OASL", "CCL3", "DDX58", "VCAM1", "PKDCC", "TDO2", "FOXF2", "NGFR", "COL9A3", "S100B", "ITGA6", "GAS7", "MLANA", "SOX10", "MITF", "PMEL", "AXL", "SERPINE1", "BGN")
write.table(fm10_barcoded.markers, file=paste0(plotDirectory,'scTransformMarkers_snn06.tsv'), col.names = TRUE, sep='\t')
plot = ggplot(data = filter(fm10_barcoded.markers, !gene %in% list), aes(y=avg_log2FC, x = cluster), color = "gray93", size = 1.5, shape = 16) +
  theme_classic() +
  geom_jitter(width = 0.2, shape = 16)+
  geom_jitter(data = filter(fm10_barcoded.markersSubset, gene %in% list), aes(x=cluster, y=avg_log2FC), color = "red", width = 0.2, shape = 16) +
  geom_text_repel(data = filter(fm10_barcoded.markersSubset, gene %in% list), aes(x=cluster, y=avg_log2FC, label = gene), color = "red") +
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plotDirectory, 'FM10_seuratAllClusterGenesWBGN.svg'), width = 12, height = 6)


#Gene painting
genes <- c("ACTA2","IFIT2", "VCAM1", "NGFR", "MITF", "AXL")
DefaultAssay(fm10_barcoded) <- "SCT"
plot <- FeaturePlot(fm10_barcoded, features = genes, slot = "data" )
ggsave(plot = plot, filename = paste0(plotDirectory,"FM10_FeaturePlot_genePainting.svg"), width = 10, height = 10)

#UMAP for clones present in both samples/twins
barcode_twins <- vector()

for (i in unique(barcodes$BC50StarcodeD8)) {
  if(all(c("S1","S2") %in% barcodes$sampleNum[barcodes$BC50StarcodeD8 ==i ] )){
    barcode_twins <- c(barcode_twins, i) 
  }
}

fm10_barcoded$Sample <- substr(colnames(fm10_barcoded),1,2)
randSampBarcodes <-barcode_twins
plots <- list()
for (i in randSampBarcodes) {
  cell_names_S1  <-colnames(fm10_barcoded[,fm10_barcoded$barcodes == i & fm10_barcoded$Sample=="S1"])
  cell_names_S2  <-colnames(fm10_barcoded[,fm10_barcoded$barcodes == i & fm10_barcoded$Sample=="S2"])
  plots[[which(randSampBarcodes == i)]] <- DimPlot(fm10_barcoded, cells.highlight = list(cell_names_S1, cell_names_S2), 
                                                   cols.highlight = c("coral2", "darkmagenta")) + ggtitle(i) +NoLegend()
} 
plot <- ggarrange(plotlist = plots, nrow = 5, ncol =5)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM10_Twin_memory.svg"), width = 20, height = 20)


#Sepcific twin clone examples

randSampBarcodes <- randSampBarcodes[c(2,5,7,6)]
plots <- list()
for (i in randSampBarcodes) {
  cell_names_S1  <-colnames(fm10_barcoded[,fm10_barcoded$barcodes == i & fm10_barcoded$Sample=="S1"])
  cell_names_S2  <-colnames(fm10_barcoded[,fm10_barcoded$barcodes == i & fm10_barcoded$Sample=="S2"])
  plots[[which(randSampBarcodes == i)]] <- DimPlot(fm10_barcoded, cells.highlight = list(cell_names_S1, cell_names_S2), 
                                                   cols.highlight = c("coral2", "darkmagenta")) + ggtitle(i) + NoLegend()
} 
plot <- ggarrange(plotlist = plots, nrow = 1)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM10_Twin_memory_specificClones.svg"), width = 12, height =3)

#select barcodes with >10 cells
fm10_rand <- fm10_barcoded[,fm10_barcoded$nSize > 10]
plots <- list()
for (bc_samp in unique(fm10_rand$barcodes)) {
  cell_names <- colnames(fm10_rand[,fm10_rand$barcodes == bc_samp])
  plots[[which(unique(fm10_rand$barcodes) == bc_samp)]] <- DimPlot(fm10, cells.highlight = cell_names  )+ NoLegend()
}
plot <- ggarrange(plotlist = plots,  ncol=5, nrow = 5)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM10_UMAP_cloneSize_10+.svg"), width = 10, height = 10)


#Barcodes for 3-10 clone siz

barcodes <- list()
for(i in 2:10){
  if(sum(fm10_barcoded$nSize == i)>1){
    fm10_rand <- fm10_barcoded[,fm10_barcoded$nSize == i]
    plots <- list()
    barcodes_size <- unique(fm10_rand$barcodes)
    n <- length(barcodes_size)
    n <- max((floor(n/5)),1)
    n <- min(n,5)
    set.seed(100)
    barcodes_size <- sample(barcodes_size, ifelse(n==1,length(barcodes_size),(n*5)))
    barcodes[[i-1]] <- barcodes_size
    for (bc_samp in barcodes_size) {
      cell_names <- colnames(fm10_rand[,fm10_rand$barcodes == bc_samp])
      plots[[which(barcodes_size == bc_samp)]] <- DimPlot(fm10, cells.highlight = cell_names  ) + NoLegend()
    }
    plot <- ggarrange(plotlist = plots,  ncol=5, nrow = ceiling(length(barcodes_size)/5)) 
    ggsave(plot = plot, filename = paste0(plotDirectory,"FM10_UMAP_cloneSize", i, ".svg"), width = 10, height = (n*2))
  }
}


#Highlighting 5 random barcodes
barcodes <- c(barcodes[[1]][c(2,8,23)], barcodes[[2]][1], barcodes[[6]][1])
plots <- list()
for (i in 1:length(barcodes)) {
  cell_names <- colnames(fm10_barcoded[,fm10_barcoded$barcodes == barcodes[i]])
  plots[[i]] <- DimPlot(fm10, cells.highlight = cell_names , cols.highlight = "hotpink3") + NoLegend() + ggtitle(barcodes[i])
  
}
plot <- ggarrange(plotlist = plots, nrow = 1)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM10_UMAP_specificClones.svg"), width =15, height =3 )
