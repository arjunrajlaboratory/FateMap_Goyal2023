########################################################################################
#Description: Clone painting of FM08
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
library(viridis)

dataDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM08/jointAnalysis/seurat/"
plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM08/jointAnalysis/seuratOutput/"

#Load data and barcodes
fm08 <- readRDS(paste0(dataDirectory, "s1s2ScTransform_50pcs_filter.rds"))
barcodes <- read.table(paste0(dataDirectory,'barcodeCellID.tsv'), stringsAsFactors=F, header = T, sep="\t")
barcodes$cellID <- paste0(barcodes$sampleNum, "_", barcodes$cellID)

#Get sizes of barcodes
barcodes <- barcodes %>%
  group_by(BC50StarcodeD8) %>%
  mutate(nSize = length(cellID))

fm08 <- RenameCells(fm08, new.names = substr(colnames(fm08),1,nchar(colnames(fm08))-2))
fm08_barcoded <- fm08[,colnames(fm08) %in% barcodes$cellID]
rownames(barcodes) <- barcodes$cellID
barcodes <- barcodes[colnames(fm08_barcoded),]
fm08_barcoded$barcodes <- barcodes$BC50StarcodeD8
fm08_barcoded$nSize <- barcodes$nSize

#UMAP for singlets
fm08_barcoded$Singlet <- fm08_barcoded$nSize ==1
fm08_barcoded$Singlet <- factor(fm08_barcoded$Singlet , levels = c(TRUE, FALSE)) 
plot <- DimPlot(fm08_barcoded, group.by = 'Singlet' )+ 
  scale_color_manual(values= c( "hotpink3", "lightgrey")) + 
  NoLegend()
ggsave(plot = plot, filename = paste0(plotDirectory ,"FM08_UMAP_singlet.svg"), width = 10, height = 10)

#UMAP for clusters
plot <- DimPlot(fm08_barcoded)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM08_UMAP_samples.svg"), width = 10, height = 10)
fm08_barcoded <- FindNeighbors(fm08_barcoded)
fm08_barcoded <- FindClusters(object=fm08_barcoded, resolution = 0.6, verbose = FALSE)
plot <- DimPlot(fm08_barcoded)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM08_UMAP_seuratCluster.svg"), width = 10, height = 10)

#Gene markers for clusters
fm08_barcoded.markers <- FindAllMarkers(fm08_barcoded, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fm08_barcoded.markersSubset <- fm08_barcoded.markers %>% filter(p_val_adj <0.05)
list = c("ACTG2", "ACTA2", "MYOCD", "TAGLN", "IFIT2", "OASL", "CCL3", "DDX58", "VCAM1", "PKDCC", "TDO2", "FOXF2", "NGFR", "COL9A3", "S100B", "ITGA6", "GAS7", "MLANA", "SOX10", "MITF", "PMEL", "AXL", "SERPINE1", "BGN")
write.table(fm08_barcoded.markers, file=paste0(plotDirectory,'scTransformMarkers_snn06.tsv'), col.names = TRUE, sep='\t')
plot = ggplot(data = filter(fm08_barcoded.markers, !gene %in% list), aes(y=avg_log2FC, x = cluster), color = "gray93", size = 1.5, shape = 16) +
  theme_classic() +
  geom_jitter(width = 0.2, shape = 16)+
  geom_jitter(data = filter(fm08_barcoded.markersSubset, gene %in% list), aes(x=cluster, y=avg_log2FC), color = "red", width = 0.2, shape = 16) +
  geom_text_repel(data = filter(fm08_barcoded.markersSubset, gene %in% list), aes(x=cluster, y=avg_log2FC, label = gene), color = "red") +
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plotDirectory, 'FM08_seuratAllClusterGenesWBGN.svg'), width = 12, height = 6)

#Plot for top markers
topDEG <- fm08_barcoded.markers$gene[order(fm08_barcoded.markers$avg_log2FC, decreasing = T)]
#topDEG <- topDEG[1:20]
topDEG <- c("TAGLN", "S100B", "TOP2A","NEAT1", "CRHBP")
plot <- FeaturePlot(fm08_barcoded, features = topDEG, slot = "data" )
ggsave(plot = plot, filename = paste0(plotDirectory,"FM08_FeaturePlot_genePainting_topDEG.svg"), width = 7, height = 10)


#Gene painting
genes <- c("ACTA2","IFIT2", "VCAM1", "NGFR", "MITF", "AXL")
DefaultAssay(fm08_barcoded) <- "SCT"
plot <- FeaturePlot(fm08_barcoded, features = genes, slot = "data" )
ggsave(plot = plot, filename = paste0(plotDirectory,"FM08_FeaturePlot_genePainting.svg"), width = 7, height = 10)

genes <- c("MITF")
DefaultAssay(fm08_barcoded) <- "SCT"
plot <- FeaturePlot(fm08_barcoded, features = genes, slot = "data" , max.cutoff = 2)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM08_FeaturePlot_genePainting_MITF_cutOff2.svg"), width = 3, height = 3)
DefaultAssay(fm08_barcoded) <- "SCT"
plot <- FeaturePlot(fm08_barcoded, features = genes, slot = "data" , max.cutoff = 1)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM08_FeaturePlot_genePainting_MITF_cutOff1.svg"), width = 3, height = 3)

#UMAP for clones present in both samples/twins
barcode_twins <- vector()

for (i in unique(barcodes$BC50StarcodeD8)) {
  if(all(c("S1","S2") %in% barcodes$sampleNum[barcodes$BC50StarcodeD8 ==i ] )){
    barcode_twins <- c(barcode_twins, i) 
  }
}

fm08_barcoded$Sample <- substr(colnames(fm08_barcoded),1,2)
set.seed(100)
randSampBarcodes <-barcode_twins
plots <- list()
for (i in randSampBarcodes) {
  cell_names_S1  <-colnames(fm08_barcoded[,fm08_barcoded$barcodes == i & fm08_barcoded$Sample=="S1"])
  cell_names_S2  <-colnames(fm08_barcoded[,fm08_barcoded$barcodes == i & fm08_barcoded$Sample=="S2"])
  plots[[which(randSampBarcodes == i)]] <- DimPlot(fm08_barcoded, cells.highlight = list(cell_names_S1, cell_names_S2), 
                                                   cols.highlight = c("coral2", "darkmagenta")) + ggtitle(i) +NoLegend()
} 
plot <- ggarrange(plotlist = plots, nrow = 5, ncol =5)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM08_Twin_memory.svg"), width = 20, height = 20)


#Sepcific twin clone examples

randSampBarcodes <- randSampBarcodes[c(3,7,17,20, 9)]
plots <- list()
for (i in randSampBarcodes) {
  cell_names_S1  <-colnames(fm08_barcoded[,fm08_barcoded$barcodes == i & fm08_barcoded$Sample=="S1"])
  cell_names_S2  <-colnames(fm08_barcoded[,fm08_barcoded$barcodes == i & fm08_barcoded$Sample=="S2"])
  plots[[which(randSampBarcodes == i)]] <- DimPlot(fm08_barcoded, cells.highlight = list(cell_names_S1, cell_names_S2), 
                                                   cols.highlight = c("coral2", "darkmagenta")) + ggtitle(i) + NoLegend()
} 
plot <- ggarrange(plotlist = plots, nrow = 1)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM08_Twin_memory_specificClones.svg"), width = 15, height =4)

#select barcodes with >10 cells
fm08_rand <- fm08_barcoded[,fm08_barcoded$nSize > 10]
plots <- list()
for (bc_samp in unique(fm08_rand$barcodes)) {
  cell_names <- colnames(fm08_rand[,fm08_rand$barcodes == bc_samp])
  plots[[which(unique(fm08_rand$barcodes) == bc_samp)]] <- DimPlot(fm08, cells.highlight = cell_names  )+ NoLegend()
}
plot <- ggarrange(plotlist = plots,  ncol=5, nrow = 5)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM08_UMAP_cloneSize_10+.svg"), width = 10, height = 10)


# #Barcodes for 3-10 clone siz
# 
# barcodes <- list()
# for(i in 3:10){
#   fm08_rand <- fm08_barcoded[,fm08_barcoded$nSize == i]
#   plots <- list()
#   barcodes_size <- unique(fm08_rand$barcodes)
#   n <- length(barcodes_size)
#   n <- max((floor(n/5)),1)
#   n <- min(n,5)
#   set.seed(100)
#   barcodes_size <- sample(barcodes_size, ifelse(n==1,length(barcodes_size),(n*5)))
#   barcodes[[i-2]] <- barcodes_size
#   for (bc_samp in barcodes_size) {
#     cell_names <- colnames(fm08_rand[,fm08_rand$barcodes == bc_samp])
#     plots[[which(barcodes_size == bc_samp)]] <- DimPlot(fm08, cells.highlight = cell_names  ) + NoLegend()
#   }
#   plot <- ggarrange(plotlist = plots,  ncol=5, nrow = n) + ggtitle(paste0("Clone size = ", i))
#   ggsave(plot = plot, filename = paste0(plotDirectory,"FM08_UMAP_cloneSize", i, ".svg"), width = 10, height = 10)
# }


#Highlighting 5 random barcodes
barcodes <- unique(fm08_barcoded$barcodes[fm08_barcoded$nSize > 10])[c(2,3,4,7,8)]
barcodes <- c(barcodes,fm08_barcoded$barcodes[fm08_barcoded$Singlet=="TRUE"][sample(1:sum(fm08_barcoded$Singlet=="TRUE"), 1)] )
plots <- list()
for (i in 1:length(barcodes)) {
  cell_names <- colnames(fm08_barcoded[,fm08_barcoded$barcodes == barcodes[i]])
  plots[[i]] <- DimPlot(fm08, cells.highlight = cell_names , cols.highlight = "hotpink3") + NoLegend() + ggtitle(barcodes[i])
  
}
plot <- ggarrange(plotlist = plots, nrow = 1)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM08_UMAP_specificClones.svg"), width = 18, height = 3)




######################
# 
# create_lpr_theme <- function(){
#   lpr_theme <- ggplot2::theme_bw() + ggplot2::theme_classic() +
#     ggplot2::theme(text = ggplot2::element_text(size = 32),
#                    plot.title = ggplot2::element_text(face = "bold",
#                                                       hjust = 0.5, size = ggplot2::rel(0.5)), axis.text = ggplot2::element_text(size = ggplot2::rel(0.8),
#                                                                                                                                 color = "black", face = "bold"), axis.title = ggplot2::element_text(size = ggplot2::rel(0.6),
#                                                                                                                                                                                                     face = "bold"), legend.title = ggplot2::element_blank(),
#                    legend.position = "bottom", axis.line = element_line(size = 2),
#                    axis.ticks = element_line(size = 2), strip.background = element_rect(size = 2),
#                    strip.text = element_text(size = ggplot2::rel(0.7),
#                                              face = "bold"))
#   return(lpr_theme)
# }
# 
# 
# #Gene painting
# logNormalizeddata <- as.data.frame(t(as.matrix(fm08_barcoded[['SCT']]@data)))
# logNormalizeddataUMAP <- cbind(logNormalizeddata , fm08_barcoded@reductions$umap@cell.embeddings[rownames(logNormalizeddata),])
# ggplot(logNormalizeddataUMAP, aes(x=UMAP_1, y= UMAP_2, color = MITF))+
#   geom_point( size = 1.5, shape = 16)+
#   create_lpr_theme() +
#   theme(legend.position = "bottom",
#         legend.title = element_text(size = rel(0.6)),
#         legend.text = element_text(size = rel(0.6), angle = 30),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())+
#   scale_color_distiller(palette = "Purples", direction = "horizontal")+
#   labs(color =  "Scaled\n Log UMI counts") 
