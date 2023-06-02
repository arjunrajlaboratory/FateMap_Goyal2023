########################################################################################
#Description: Clone painting of FM07
#Author: Maalavika Pillai
#Version: 1
#Edited on: 21st November 2022
########################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(svglite) 
library(ggrepel)
library(viridis)
library(VennDiagram)

dataDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM07/jointAnalysis/seurat/25_15_mt_"
data1Directory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM07/jointAnalysis/seurat/"
plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM07/jointAnalysis/seuratOutput/25_15_mt_"

#Load data and barcodes
fm07 <- readRDS(paste0(dataDirectory, "s1s2ScTransform_50pcs_filter.rds"))
barcodes <- read.table(paste0(data1Directory,'barcodeCellID2.tsv'), stringsAsFactors=F, header = T, sep="\t")
barcodes$cellID <- paste0(barcodes$sampleNum, "_", barcodes$cellID)

#Get sizes of barcodes
fm07 <- RenameCells(fm07, new.names = substr(colnames(fm07),1,nchar(colnames(fm07))-2))
fm07_barcoded <- fm07[,colnames(fm07) %in% barcodes$cellID]
rownames(barcodes) <- barcodes$cellID
barcodes <- barcodes[colnames(fm07_barcoded),]

barcodes <- barcodes %>%
  group_by(BC50StarcodeD8) %>%
  mutate(nSize = length(cellID))

fm07_barcoded$barcodes <- barcodes$BC50StarcodeD8
fm07_barcoded$nSize <- barcodes$nSize
fm07_barcoded$Sample <- substr(colnames(fm07_barcoded),1,2)
n_barcodes <- length(unique(barcodes$BC50StarcodeD8))
n_barcodes_S1 <- length(unique(barcodes$BC50StarcodeD8[barcodes$sampleNum=="S1"]))
n_barcodes_S2 <-length(unique(barcodes$BC50StarcodeD8[barcodes$sampleNum=="S2"]))
fracAll_S1 <- n_barcodes_S1/n_barcodes
fracAll_S2 <- n_barcodes_S2/n_barcodes


#UMAP for singlets
fm07_barcoded$Singlet <- fm07_barcoded$nSize ==1
fm07_barcoded$Singlet <- factor(fm07_barcoded$Singlet , levels = c(TRUE, FALSE)) 
plot <- DimPlot(fm07_barcoded, group.by = 'Singlet' )+ 
  scale_color_manual(values= c( "hotpink3", "lightgrey")) + 
  NoLegend()
ggsave(plot = plot, filename = paste0(plotDirectory ,"FM07_UMAP_singlet.svg"), width = 10, height = 10)

#UMAP for clusters
plot <- DimPlot(fm07_barcoded)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM07_UMAP_samples.svg"), width = 10, height = 10)
svglite(paste0(plotDirectory,"SampleDistribution_Pie.svg"))
pie(c(sum(fm07_barcoded$Sample=="S1"),sum(fm07_barcoded$Sample=="S2")), 
    labels = c("Continuous", "Discontinuous"), col = c("#EFA282", "#87DED8"))
dev.off()

fm07_barcoded <- FindNeighbors(fm07_barcoded)
fm07_barcoded <- FindClusters(object=fm07_barcoded, resolution = 0.6, verbose = FALSE)
plot <- DimPlot(fm07_barcoded)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM07_UMAP_seuratCluster.svg"), width = 10, height = 10)

#Gene markers for clusters
fm07_barcoded.markers <- FindAllMarkers(fm07_barcoded, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fm07_barcoded.markersSubset <- fm07_barcoded.markers %>% filter(p_val_adj <0.05)
list = c("ACTG2", "ACTA2", "MYOCD", "TAGLN", "IFIT2", "OASL", "CCL3", "DDX58", "VCAM1", "PKDCC", "TDO2", "FOXF2", "NGFR", "COL9A3", "S100B", "ITGA6", "GAS7", "MLANA", "SOX10", "MITF", "PMEL", "AXL", "SERPINE1", "BGN")
write.table(fm07_barcoded.markers, file=paste0(plotDirectory,'scTransformMarkers_snn06.tsv'), col.names = TRUE, sep='\t')
plot = ggplot(data = filter(fm07_barcoded.markers, !gene %in% list), aes(y=avg_log2FC, x = cluster), color = "gray93", size = 1.5, shape = 16) +
  theme_classic() +
  geom_jitter(width = 0.2, shape = 16)+
  geom_jitter(data = filter(fm07_barcoded.markersSubset, gene %in% list), aes(x=cluster, y=avg_log2FC), color = "red", width = 0.2, shape = 16) +
  geom_text_repel(data = filter(fm07_barcoded.markersSubset, gene %in% list), aes(x=cluster, y=avg_log2FC, label = gene), color = "red") +
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plotDirectory, 'FM07_seuratAllClusterGenesWBGN.svg'), width = 12, height = 6)

#Plot for top markers
topDEG <- fm07_barcoded.markers$gene[order(fm07_barcoded.markers$avg_log2FC, decreasing = T)]
#topDEG <- topDEG[1:20]
topDEG <- c("TAGLN", "S100B", "TOP2A","NEAT1", "CRHBP")
plot <- FeaturePlot(fm07_barcoded, features = topDEG, slot = "data" )
ggsave(plot = plot, filename = paste0(plotDirectory,"FM07_FeaturePlot_genePainting_topDEG.svg"), width = 7, height = 10)


#Gene painting
genes <- c("ACTA2","IFIT2", "VCAM1", "NGFR", "MLANA", "AXL")
DefaultAssay(fm07_barcoded) <- "SCT"
plot <- FeaturePlot(fm07_barcoded, features = genes, slot = "data" )
ggsave(plot = plot, filename = paste0(plotDirectory,"FM07_FeaturePlot_genePainting.svg"), width = 7, height = 10)


#Histogram for MLANA
logNormalizeddata <- as.data.frame(t(as.matrix(fm07_barcoded[['SCT']]@data)))
hist(logNormalizeddata$MLANA, breaks = 10)
fm07_barcoded$MLANA_cluster <- logNormalizeddata$MLANA>2
MLANA_cluster_data <- data.frame("MLANA positive"= fm07_barcoded$MLANA_cluster, "Singlet" = fm07_barcoded$Singlet,
                                 "Sample" = fm07_barcoded$Sample, "nSize" = fm07_barcoded$nSize , 
                                 "barcode" = fm07_barcoded$barcodes)
MLANA_cluster_data <- MLANA_cluster_data %>% 
  filter(MLANA.positive) %>%
  group_by(barcode) %>% 
  mutate("FractionCloneMLANApos" = length(Sample)/ unique(nSize))

MLANA_cluster_data <- MLANA_cluster_data %>% 
  filter(FractionCloneMLANApos >0.5) 

MLANA_cluster_singlet <- MLANA_cluster_data %>%
  group_by(Sample, barcode) %>%
  summarize("Singlet" = sum(as.logical(Singlet)), "nSize" = nSize)
 
MLANA_cluster_singlet <- MLANA_cluster_singlet %>%
  group_by(Sample) %>%
  summarize("Singlet" = sum(as.logical(Singlet))/length(Singlet), avgSize = mean(nSize)) %>%
  as.data.frame()

svglite(paste0(plotDirectory, "MLANApos_SingletonFraction_pie.svg"), width = 10, height = 6.8)
par(mfrow = c(1,2))
pie(c(MLANA_cluster_singlet[1,2], 1- MLANA_cluster_singlet[1,2]), labels =  c("Singletons", "Larger colonies"), main = "Continuous treatment")
pie(c(MLANA_cluster_singlet[2,2], 1- MLANA_cluster_singlet[2,2]), labels =  c("Singletons", "Larger colonies"), main = "Discontinuous treatment")
dev.off()

#UMAP for clones present in both samples/twins
barcode_twins <- vector()

for (i in unique(barcodes$BC50StarcodeD8)) {
  if(all(c("S1","S2") %in% barcodes$sampleNum[barcodes$BC50StarcodeD8 ==i ] )){
    barcode_twins <- c(barcode_twins, i) 
  }
}
n_overlap <- length(barcode_twins)

#Plot for overlaping barcodes
cell_names_S1  <-colnames(fm07_barcoded[,fm07_barcoded$barcodes %in% barcode_twins & fm07_barcoded$Sample=="S1"])
cell_names_S2  <-colnames(fm07_barcoded[,fm07_barcoded$barcodes%in% barcode_twins & fm07_barcoded$Sample=="S2"])
plot <- DimPlot(fm07_barcoded, cells.highlight = list(cell_names_S1, cell_names_S2), 
                  cols.highlight = c("coral2", "darkmagenta")) + NoLegend()+
  annotate("text", x=8, y=5, label= paste0("Total number of barcodes = ", n_barcodes, size = 6)) + 
  annotate("text", x=8, y=4.7, label= paste0("Fraction of all barcodes in continuous treatment = ", round(fracAll_S1,2), size=6)) +
  annotate("text", x=8, y=4.4, label= paste0("Fraction of all barcodes in discontinuous treatment = ", round(fracAll_S2,2), size=6)) +
  annotate("text", x=8, y=4.1, label= paste0("Fraction of continuous barcodes overlapping with discontinuous = ", round((n_overlap/n_barcodes_S1),2), size=6)) +
  annotate("text", x=8, y=3.8, label= paste0("Fraction of discontinuous barcodes overlapping with continuous = ", round((n_overlap/n_barcodes_S2),2), size=6)) +
  xlim(c(-5,12))
ggsave(plot = plot, filename = paste0(plotDirectory,"FM07_OverlappingBarcodes.svg"), width = 12.5, height = 10)

#Piechart for fractions
n_barcodes_ind <- c(n_barcodes_S1-n_overlap, n_barcodes_S2-n_overlap, n_overlap)
svglite(paste0(plotDirectory,"BarcodeOverlap_Venn.svg"))
draw.pairwise.venn(n_barcodes_S1, n_barcodes_S2, n_overlap,
                   category=c("Continuous","Discontinuos"),
                   fill=c("darkmagenta","coral2"))
dev.off()


#Twin memory for all clones
set.seed(100)
randSampBarcodes <- sample(barcode_twins,25)
#132,107,120,81,85,91,55,70,

plots <- list()
for (i in randSampBarcodes) {
  cell_names_S1  <-colnames(fm07_barcoded[,fm07_barcoded$barcodes == i & fm07_barcoded$Sample=="S1"])
  cell_names_S2  <-colnames(fm07_barcoded[,fm07_barcoded$barcodes == i & fm07_barcoded$Sample=="S2"])
  plots[[which(randSampBarcodes == i)]] <- DimPlot(fm07_barcoded, cells.highlight = list(cell_names_S1, cell_names_S2), 
                                                   cols.highlight = c("coral2", "darkmagenta")) + ggtitle(i) +NoLegend()
} 
plot <- ggarrange(plotlist = plots, nrow = 5, ncol =5)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM07_Twin_memory.svg"), width = 20, height = 20)

#Sepcific twin clone examples
#randSampBarcodes <- barcode_twins[c(130,128, 127, 112, 109, 81, 114, 117,123,108,45,43 )]
randSampBarcodes <- barcode_twins[c(130, 112, 114, 108,109, 128, 45,43)]
plots <- list()
for (i in randSampBarcodes) {
  cell_names_S1  <-colnames(fm07_barcoded[,fm07_barcoded$barcodes == i & fm07_barcoded$Sample=="S1"])
  cell_names_S2  <-colnames(fm07_barcoded[,fm07_barcoded$barcodes == i & fm07_barcoded$Sample=="S2"])
  plots[[which(randSampBarcodes == i)]] <- DimPlot(fm07_barcoded, cells.highlight = list(cell_names_S1, cell_names_S2), 
                                                   cols.highlight = c("coral2", "darkmagenta")) + ggtitle(i) +NoLegend()
} 
plot <- ggarrange(plotlist = plots, nrow = 5, ncol =5)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM07_Twin_specific_clones.svg"), width = 20, height = 20)

#Discontinuous clones specific cases
# randSampBarcodes <- c("ATTGATCCAGTTCAACCAGTTCGTCAAGCAGGTCGTCTTCCACTACCTCT","ATTGTTGCTGTTGGTGTTCGAGGTCCAGCACTACTTGATGTACAACGAGT",
#                       "ATACTTGTACGTGAACTAGGTGCAGCAGCTGTACTACCTGTAGTTGTTCG","AGAGGTGTTGGAGGTCTTCTACTAGATGTAATTGGTCTTCATGTACCTCG",
#                       "ATTGTTGATCATCCAGTTAATGTTGGTGTTGGAGATCTTGATCCTGCTGG","ATAGGACGTGCTGCTGTTCATCGTGAAGGTGGTGTTCGTCTTGGTGGTCA",
#                       "ATTCTTCCTCCTCTTCGACGAGCAGCTCCTCATGCTCCTGTACTACGAGT","AGTCAAGTAGTAGTTGGAGTTGCTGAAGGACGAGAAGGTGAAGCTGGAGA",
#                       "ATTGCTGGTGATGGAGATGGTCAACCAGATCTTCAACCAGCTCCACTAGA","ATACCTCATGGAGGTCTTCTAGATGAAGAACCTGATGTTCTTGAAGGACA",
#                       "ATTGGAGGTCTTCTACTACAAGCTCTAGGTCCTGAACCTGCAGTACTAGA","ATTGGTGCTCATGGACCTGAACGTCATCCTCCAGCTCATCTACGTCCTGC")
randSampBarcodes <- c("ATAGGACGTGCTGCTGTTCATCGTGAAGGTGGTGTTCGTCTTGGTGGTCA","ATACCTCATGGAGGTCTTCTAGATGAAGAACCTGATGTTCTTGAAGGACA",
                      "ATTGATCCAGTTCAACCAGTTCGTCAAGCAGGTCGTCTTCCACTACCTCT","ATTCTTCCTCCTCTTCGACGAGCAGCTCCTCATGCTCCTGTACTACGAGT",
                      "AGAGGTGTTGGAGGTCTTCTACTAGATGTAATTGGTCTTCATGTACCTCG", "ATTGGAGGTCTTCTACTACAAGCTCTAGGTCCTGAACCTGCAGTACTAGA",
                      "ATTGGTGCTCATGGACCTGAACGTCATCCTCCAGCTCATCTACGTCCTGC")
plots <- list()
for (i in randSampBarcodes) {
  cell_names<-colnames(fm07_barcoded[,fm07_barcoded$barcodes == i & fm07_barcoded$Sample=="S2"])
  plots[[which(randSampBarcodes == i)]] <- DimPlot(fm07_barcoded, cells.highlight = list(cell_names), 
                                                   cols.highlight = c("hotpink3")) + ggtitle(i) +NoLegend()
} 
plot <- ggarrange(plotlist = plots, nrow = 5, ncol =5)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM07_Discontinuous_specific_clones.svg"), width = 20, height = 20)


##########
#Total unique barcodes for each sample
###Temp for FateMap
continuous = 2510
discontinuous = 703
a = as_tibble(c(continuous,discontinuous)) %>% mutate(type = c("discontinuous", "continuous" ))

plot_FateMap = ggplot (a, aes(type, value))+
  geom_point() +
  geom_path(group =1) +
  ylim(0,2600) +
  theme_classic()

plotDirectory = "/Users/ygy1258/Dropbox (RajLab)/FateMap/REVISIONS/analysisPlots/FM07/"

ggsave(plot_FateMap, file = paste0(plotDirectory, 'totalUniqueBarcodes.svg'), width = 2, height = 4)
