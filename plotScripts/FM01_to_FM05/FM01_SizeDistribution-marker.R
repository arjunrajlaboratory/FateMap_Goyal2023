#Description: ANOVA F-test to identify DEGs in clones
#Aim:Create  histogram for clone size, identify DEGs in top n clones, compare to DEGs
#from clusters
#Edited to add plot without singlets
#Author: Maalavika Pillai
#Version: 2
#Created on: 12/22/22
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
library(corrplot)
library(rstatix)
library(ggpubr)

# Load the relevant datasets (filteredbarcode matrix files)

homeDirectory <-"~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM01/Data/"
plotDirectory <-"~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM01/Plots/"

##FM01_sct_integratedData_sample1sample2
fm01<- readRDS(file = paste0(homeDirectory,'s1s2ScTransform_50pcs_filter.rds'))  
barcodes <- read.table(paste0(homeDirectory,'barcodeCellID.tsv'), stringsAsFactors=F, header = T, sep="\t")
#for S2 read this barcodes file
#barcodes <- read.table(paste0(homeDirectory,'barcodeCellID_S2.tsv'), stringsAsFactors=F, header = T, sep="\t")
barcodes$cellID <- paste0(barcodes$sampleNum, "_", barcodes$cellID)

#Get sizes of barcodes
barcodes <- barcodes %>%
  group_by(BC50StarcodeD8) %>%
  mutate(nSize = length(cellID))

fm01_barcoded <- fm01[,colnames(fm01) %in% barcodes$cellID]
rownames(barcodes) <- barcodes$cellID
barcodes <- barcodes[colnames(fm01_barcoded),]
fm01_barcoded$barcodes <- barcodes$BC50StarcodeD8
fm01_barcoded$nSize <- barcodes$nSize
fm01_barcoded$CloneSize <- ifelse(fm01_barcoded$nSize==1, "Singleton", "Small Colony")
fm01_barcoded$CloneSize[fm01_barcoded$nSize>3] <- "Large Colony"
barcodes$CloneSize <- fm01_barcoded$CloneSize


#Plot distribution for various clusters
barcodes$cluster <- fm01_barcoded$integrated_snn_res.0.6

#Test run
# barcodes_sub <- barcodes %>% filter(cluster %in% c(7,11,8))
# ggplot(barcodes_sub, aes(x = cluster, fill = CloneSize)) + 
#   geom_bar(position = "fill")

#We use chi square test to see if the cluster is independent of colony size or not
barcodes_sub <- barcodes %>% filter(!cluster %in% c(5,9,10,12:14))
barcodes_sub$cluster <- as.character(barcodes_sub$cluster)
barcodes_sub$cluster[barcodes_sub$cluster %in% c(0,3)] <- "MLANA"
barcodes_sub$cluster[barcodes_sub$cluster %in% c(7)] <- "NGFR"
#barcodes_sub$cluster[barcodes_sub$cluster %in% c(11)] <- "IFIT2"
barcodes_sub$cluster[barcodes_sub$cluster %in% c(8)] <- "ACTA2"
barcodes_sub$cluster[barcodes_sub$cluster %in% c(1,2,11)] <- "AXL"
barcodes_sub$cluster[barcodes_sub$cluster %in% c(6,4,15)] <- "VCAM1"

barcodes_sub$cluster <- factor(barcodes_sub$cluster, levels = c("MLANA", "NGFR", "AXL", "VCAM1","ACTA2"))

barcodes_sub_colonies <- barcodes_sub %>%
  group_by(BC50StarcodeD8, cluster) %>%
  summarize(nSize = length(BC50StarcodeD8))

#t-test to check if difference in colony size is statistically significant
stat.test <- barcodes_sub_colonies %>%
  ungroup() %>%
  wilcox_test(nSize~cluster)
stat.test <- stat.test %>% add_xy_position(x = "cluster")
stat.test <- stat.test[order(stat.test$y.position, decreasing = F), ]
stat.test <- stat.test %>%
  mutate(y.position = seq(20,50,length.out=length(y.position)))
#Plot average size of colonies for each marker
set.seed(100)
plot <- ggplot(barcodes_sub_colonies, aes(x = cluster, y = nSize, fill = cluster))+
  geom_boxplot(width = 0.5, alpha =0.5)+
  stat_pvalue_manual(stat.test, tip.length = 0.02,bracket.nudge.y = 0, size = 8,inherit.aes = FALSE)+
  theme_classic((base_size = 28))+
  ylab("Size of colonies") + 
  xlab("Marker")+ 
  ylim(c(0,50))
ggsave(plot = plot, filename = paste0(plotDirectory, "SizeDistribution_markers_axis_cut.svg"), width = 10, height = 10)

stat.test <- barcodes_sub_colonies %>%
  ungroup() %>%
  wilcox_test(nSize~cluster)
stat.test <- stat.test %>% add_xy_position(x = "cluster")
stat.test <- stat.test[order(stat.test$y.position, decreasing = F), ]
stat.test <- stat.test %>%
  mutate(y.position = seq(30,100,length.out=length(y.position)))


set.seed(100)
plot <- ggplot(barcodes_sub_colonies, aes(x = cluster, y = nSize, fill = cluster))+
  geom_boxplot(width = 0.5, alpha =0.5)+
  stat_pvalue_manual(stat.test, tip.length = 0.02,bracket.nudge.y = 0, size = 8,inherit.aes = FALSE)+
  theme_classic((base_size = 28))+
  ylab("Size of colonies") + 
  xlab("Marker")
ggsave(plot = plot, filename = paste0(plotDirectory, "SizeDistribution_markers_axis_notcut.svg"), width = 10, height = 10)

barcodes_sub_colonies <- barcodes_sub_colonies %>%
  ungroup() %>%
  filter(!nSize==1)
stat.test <- barcodes_sub_colonies %>%
  ungroup() %>%
  wilcox_test(nSize~cluster)
stat.test <- stat.test %>% add_xy_position(x = "cluster")
stat.test <- stat.test[order(stat.test$y.position, decreasing = F), ]
stat.test <- stat.test %>%
  mutate(y.position = seq(30,100,length.out=length(y.position)))

set.seed(100)
plot <- ggplot(barcodes_sub_colonies, aes(x = cluster, y = nSize, fill = cluster))+
  geom_boxplot(width = 0.5, alpha =0.5)+
  stat_pvalue_manual(stat.test, tip.length = 0.02,bracket.nudge.y = 0, size = 8,inherit.aes = FALSE)+
  theme_classic((base_size = 28))+
  ylab("Size of colonies") + 
  xlab("Marker")
ggsave(plot = plot, filename = paste0(plotDirectory, "SizeDistribution_markers_noSingleton_axis_notcut.svg"), width = 10, height = 10)

#For S2 
#ggsave(plot = plot, filename = paste0(plotDirectory, "SizeDistribution_markers_S2.svg"), width = 10, height = 10)

barcodes_sub <- barcodes_sub %>%
  ungroup() %>%
  select(CloneSize, cluster)%>%
  group_by(cluster,CloneSize) %>%
  mutate("ColonySize" = length(CloneSize)) 
barcodes_sub$ColonySize <- as.numeric(barcodes_sub$ColonySize)

barcode_sub_wide <- barcodes_sub %>%
  ungroup() %>%
  pivot_wider(names_from = CloneSize, values_from = ColonySize, values_fn = unique)
barcode_sub_wide <- as.data.frame(barcode_sub_wide)
rownames(barcode_sub_wide) <-  barcode_sub_wide$cluster
barcode_sub_wide[is.na(barcode_sub_wide)] <- 0
barcode_sub_wide <- barcode_sub_wide[,-1]

chiTest <- chisq.test(barcode_sub_wide)
svglite(filename =  paste0(plotDirectory, "PearsonResiduals_SizeDistributionMarker.svg"), width = 6, height = 10)
corrplot(chiTest$residuals, is.cor = FALSE)
dev.off()

