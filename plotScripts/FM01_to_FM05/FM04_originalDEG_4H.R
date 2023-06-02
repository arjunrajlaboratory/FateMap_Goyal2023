########################################################################################
#Description: FM04 - add statistical test to fig 4H and identify DEGs overlapping
#with cell cycle
#Author: Maalavika Pillai
#Version: 1
#Edited on: 11/2/22
########################################################################################
########################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(svglite) 
library(ggrepel)
library(VennDiagram)
library(RColorBrewer)

dataDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM04/seurat/"
plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM04/seuratOutput/"

fractionFinalClusters <- read.table(paste0(dataDirectory,"FM04_fractionFinalClusters_S1_snn05.tsv" ),sep ="\t")
fractionFinalClusters$type <- factor(fractionFinalClusters$type, levels = c("random", "colonies"))
plot = ggplot(fractionFinalClusters, aes(x = type, y=nfractionFinal, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  geom_jitter(width = 0.1, shape = 16) +
  stat_compare_means( size =4, label.x.npc = "left")+
  scale_fill_manual(values=c("grey","goldenrod3")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none") +
  ylim(0,1)
fractionFinalClusters %>% group_by(type) %>% summarize(mean(nfractionFinal))
ggsave(plot, file = paste0(plotDirectory, 'FM04_ColonySpecificPlotColony_V1_withDotsFinal_S1_05.svg'), width = 4, height = 8)

fm04 <- readRDS(paste0(dataDirectory, "s1s2ScTransform_50pcs_filter.rds"))
ccGenes <- unlist(cc.genes)
fm04<- FindNeighbors(fm04)
fm04 <- FindClusters(object=fm04, resolution = 0.5, verbose = FALSE)
fm04.markers <- FindAllMarkers(fm04)
write.table(fm04.markers,paste0(dataDirectory,"FM04_DEG_seurat clusters.tsv"), quote = F, sep = "\t")
DEG <- unique(fm04.markers$gene)
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(DEG, ccGenes),
  category.names = c( "Clonal DEGs", "Cell cycle genes" ),
  filename = paste0(plotDirectory, "CloneGenes_CellCycleGenes_Venn.png" ),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
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
