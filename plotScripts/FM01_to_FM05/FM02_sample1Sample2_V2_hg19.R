#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2020_FM02/jointAnalysis/seurat/hg19/sample1Sample2_V1_hg19.R
##############################################

library(Seurat)
library(ggplot2)
library(tidyverse)
library(fuzzyjoin)
library(pheatmap)
library(RColorBrewer)
options(future.globals.maxSize = 4000 * 1024^2)

#########
# workflow adapted from https://satijalab.org/seurat/v3.0/integration.html
#########
home1Directory <-"/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM02/"
home2Directory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM02/'
home3Directory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/gDNA_BarcodeData/FM02/dataFor10X/'
plot3Directory <- "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/plotData/analysisPlots/"
create_lpr_theme <- function(){
  lpr_theme <- ggplot2::theme_bw() + ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 32),
                   plot.title = ggplot2::element_text(face = "bold",
                                                      hjust = 0.5, size = ggplot2::rel(0.5)), axis.text = ggplot2::element_text(size = ggplot2::rel(0.8),
                                                                                                                                color = "black", face = "bold"), axis.title = ggplot2::element_text(size = ggplot2::rel(0.6),
                                                                                                                                                                                                    face = "bold"), legend.title = ggplot2::element_blank(),
                   legend.position = "bottom", axis.line = element_line(size = 2),
                   axis.ticks = element_line(size = 2), strip.background = element_rect(size = 2),
                   strip.text = element_text(size = ggplot2::rel(0.7),
                                             face = "bold"))
  return(lpr_theme)
}

s1s2Scanorama <- readRDS(paste0(home1Directory,"s1s2Scanorama_50pcs_filter_hg19.rds"))
logNormalizedCounts = as_tibble(read.table(file = paste0(home1Directory, "logNormalizedCounts_s1s2Scanorama_50pcs_filterRound_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinates = as_tibble(read.table(file = paste0(home1Directory, "umapCoordinates_s1s2Scanorama_50pcs_filter_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapClusters = as_tibble(read.table(file = paste0(home1Directory, "umapClusters_s1s2Scanorama_50pcs_filter_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))

#rm(logNormalizedCounts)
################3Figure 3E,F################################################
#####Just to plot the three with the right dimensions
s1s2Scanorama <- FindClusters(object=s1s2Scanorama,graph.name = "scanorama_snn", resolution=0.8)
DimPlot(s1s2Scanorama, reduction = "umap_scanorama", group.by = "orig.ident") + NoLegend()
clusterUMAP_08 = DimPlot(s1s2Scanorama, reduction = "umap_scanorama", label = FALSE, group.by = "scanorama_snn_res.0.8") + NoLegend()
clusterUMAP_08WCluster = DimPlot(s1s2Scanorama, reduction = "umap_scanorama", label = TRUE, label.size = 8, group.by = "scanorama_snn_res.0.8") + NoLegend()
sampleUMAP = DimPlot(s1s2Scanorama, reduction = "umap_scanorama", group.by = "orig.ident") + NoLegend()
sampleSplitUMAP_08 = DimPlot(s1s2Scanorama, reduction = "umap_scanorama",slot="RNA.data", group.by = "scanorama_snn_res.0.8", split.by = "orig.ident", label = FALSE) + NoLegend()

ggsave(clusterUMAP_08, file = paste0(plot3Directory, 'FM02_s1s2_clusterUMAP_onlyScanorma08_hg19.svg'), width = 6, height = 5.321)
ggsave(clusterUMAP_08WCluster, file = paste0(plot3Directory, 'FM02_s1s2_clusterUMAPWNUMBERS_onlyScanorma08_hg19.svg'), width = 6, height = 5.321)
ggsave(sampleUMAP, file = paste0(plot3Directory, 'FM02_s1s2_sampleUMAP_onlyScanorma_hg19.svg'), width = 6, height = 5.321)
ggsave(sampleSplitUMAP_08, file = paste0(plot3Directory, 'FM02_s1s2_sampleSplitUMAP_onlyScanorma08_hg19.svg'), width = 12, height = 5.822)

rm(s1s2Scanorama)
####here it can be done for all of markers together

plot <- ggplot(umapCoordinates, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes_string(color = logNormalizedCounts$MLANA), size = 1, shape = 16) +
  scale_color_gradient(low = "gray93", high = "darkblue") +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Scaled\n Log UMI counts")
ggsave(plot, file = paste0(plot3Directory, 'FM02_MLANA.svg'), width = 6, height = 6.429)

#############
#############################################################################
##############Barcode entries for Sample 1 and Sample 2#####################
#############################################################################
barcode50 = as_tibble(read.table(paste0(home2Directory, 'stepThreeStarcodeShavedReads.txt'), stringsAsFactors=F, header = F))
barcode50 = barcode50  %>% 
  dplyr::rename(cellID = V1,
                UMI = V2,
                BC50StarcodeD8 = V4,
                sampleNum = V8) %>%
  select(-c(V3,V5,V6,V7)) %>% unique() %>% select(-UMI) %>% filter(sampleNum %in% c("1","2"))

umiCut = 15 # Minimum UMI cutoff for reliable analysis. For the past runs I used between 3-5 for ~500k reads/sample. Here we have ~1.5mil/sample
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.  

upToLineageCounting = barcode50 %>% 
  group_by(cellID, BC50StarcodeD8, sampleNum) %>%
  summarise(nUMI = length(sampleNum)) %>%
  filter(nUMI >= umiCut) %>% 
  group_by(cellID) %>%
  mutate(nLineages = length(cellID)) 

upToLineageCounting$sampleNum = sub("^", "S", upToLineageCounting$sampleNum) 
#####
###taking only single cellid-->barcode mappings
linCountTooverlaps = upToLineageCounting %>%
  ungroup() %>%
  filter(nLineages <= linCut) %>%
  unique()

rm(barcode50)

jointUMAP = inner_join(linCountTooverlaps, umapCoordinates, by = c("cellID", "sampleNum")) %>% select(-nLineages)

#####to get cluster numbers
clusterUMAP = inner_join(umapClusters,umapCoordinates, by = c("cellID", "sampleNum"))
ggplot() +
  geom_point(data = clusterMLANAUMAP, aes(x = UMAP_1, y = UMAP_2, color = factor(scanorama_snn_res.0.8)), size = 1) +
  create_lpr_theme()

#############################################################################
##############Analysis for MLANA Cluster High Vs Low ########################
#############################################################################

clusterMLANACellIDs = umapClusters %>% filter(scanorama_snn_res.0.8 %in% c(6,5) )
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,scanorama_snn_res.0.8)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

#####basic plots of cluster 2 and Sample 1 vs Sample 2 
#####colonyBreakup sample 1 vs 2

summaryMLANAClusters = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                             colony = sum(total[nColony.y >1]))


summaryMLANAClusters = summaryMLANAClusters %>% mutate(percentColony = 100*colony/(singlet+colony))  ###used for the paper numbers

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1, shape = 16) +
  geom_point(data = clusterMLANAUMAP, aes(x = UMAP_1, y = UMAP_2, color = factor(scanorama_snn_res.0.8)), size = 1, shape = 16) +
  create_lpr_theme()+ theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM02_ClusterMLANA_S1S2.svg'), width = 6, height = 5.321)


colonyS2 = fractionColony %>% filter(nColony.y >2, sampleNum == "S2") 
colonyS2UMAP = jointUMAP %>% filter(sampleNum =="S2") %>% filter(BC50StarcodeD8 %in% c(colonyS2$BC50StarcodeD8))

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93",size = 2, shape = 16) +
  geom_point(data = colonyS2UMAP, aes(x = UMAP_1, y = UMAP_2, col=BC50StarcodeD8),size = 3, shape = 16) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM02_largeColoniesS2MLANA.svg'), width = 6, height = 5.321)

SingletS1 = fractionColony %>% filter(nColony.y ==1, sampleNum == "S1") 
selectS1UMAP = jointUMAP %>% filter(sampleNum =="S1") %>% filter(BC50StarcodeD8 %in% c(SingletS1$BC50StarcodeD8))

SingletS1All = colonySizeAll %>% filter(nColony ==1, sampleNum == "S1") ####just to check the distributions
selectS1UMAPAll = jointUMAP %>% filter(sampleNum =="S1") %>% filter(BC50StarcodeD8 %in% c(SingletS1All$BC50StarcodeD8)) ####just to check the distributions

SingletS2 = fractionColony %>% filter(nColony.y ==1, sampleNum == "S2") 
selectS2UMAP = jointUMAP %>% filter(sampleNum =="S2") %>% filter(BC50StarcodeD8 %in% c(SingletS2$BC50StarcodeD8))

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93",size = 2, shape = 16) +
  geom_point(data = selectS1UMAP, aes(x = UMAP_1, y = UMAP_2), col="hotpink3", size = 3, shape = 16) +
  geom_point(data = selectS2UMAP, aes(x = UMAP_1, y = UMAP_2), col="turquoise3",size = 3, shape = 16) +
  create_lpr_theme()
ggsave(plot, file = paste0(plot3Directory, 'FM02_SingletsS1S2_ClusterMLANA.svg'), width = 6, height = 5.321)


#############################################################################
####### singlet ONLY Cluster MLANA ######################
#############################################################################
clusterMLANACellIDs = umapClusters %>% filter(scanorama_snn_res.0.8 %in% c(6,5) )
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,scanorama_snn_res.0.8)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)

jointUMAPS1 = jointUMAP %>% filter(sampleNum=="S1")
jointUMAPS2 = jointUMAP %>% filter(sampleNum=="S2")

jointUMAPS1Barcodes = jointUMAP %>% filter(sampleNum=="S1") %>% select(BC50StarcodeD8) %>% unique()
jointUMAPS2Barcodes = jointUMAP %>% filter(sampleNum=="S2")  %>% select(BC50StarcodeD8) %>% unique()

jointBarcodesOnlyBoth = inner_join(jointUMAPS1Barcodes,jointUMAPS2Barcodes, by = "BC50StarcodeD8") 
jointBarcodesOnlyBoth = jointBarcodesOnlyBoth %>% mutate(barcodeName = c(1:nrow(jointBarcodesOnlyBoth)))
jointBarcodesOnlyBoth$barcodeName = sub("^", "B", jointBarcodesOnlyBoth$barcodeName)

jointBarcodesOnlyBothList = jointBarcodesOnlyBoth$BC50StarcodeD8

jointUMAPOnlyBoth = jointUMAP %>% filter(BC50StarcodeD8 %in% jointBarcodesOnlyBothList)
finalUMAPJoint = inner_join(jointUMAPOnlyBoth,jointBarcodesOnlyBoth, by = "BC50StarcodeD8") %>% select(-cellID, -nUMI)
fractionColonySinglet = fractionColony %>% filter(nColony.y ==1)
finalUMAPJointMLANA = inner_join(finalUMAPJoint, fractionColony, by = c("BC50StarcodeD8","sampleNum")) %>% select(-nColony.y)

finalUMAPS1 = finalUMAPJointMLANA %>% filter(sampleNum=="S1")
finalUMAPS2 = finalUMAPJointMLANA %>% filter(sampleNum=="S2")
finalUMAPS2All = finalUMAPJoint %>% filter(sampleNum=="S2")
finalUMAPS1Singlet = finalUMAPS1 %>% group_by(barcodeName,sampleNum) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony ==1)
finalUMAPS1SingletOthers = finalUMAPS1 %>% group_by(barcodeName,sampleNum) %>% summarise(nColony = length(barcodeName))
finalUMAPS2SingletSibling = finalUMAPS2 %>% group_by(barcodeName,sampleNum) %>% summarise(nColony = length(barcodeName)) 
finalUMAPS2SingletSiblingAll = finalUMAPS2All %>% group_by(barcodeName,sampleNum) %>% summarise(nColony = length(barcodeName)) 

#finalUMAPS2NonSingletSibling = finalUMAPS2 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony >1) %>% select(-nColony)

commonS1S2Singlet = inner_join(finalUMAPS1Singlet, finalUMAPS2SingletSibling, by = "barcodeName")
commonS1S2Singlet = commonS1S2Singlet$barcodeName

commonS1S2SingletAll = inner_join(finalUMAPS1Singlet, finalUMAPS2SingletSiblingAll, by = "barcodeName")
commonS1S2SingletAll = commonS1S2SingletAll$barcodeName

commonS1S2SingletNonSinglet = inner_join(finalUMAPS1SingletOthers, finalUMAPS2SingletSibling, by = "barcodeName")
commonS1S2NonSingletSinglet = commonS1S2SingletNonSinglet$barcodeName

finalUMAPJointSinglet = finalUMAPJoint %>% filter(barcodeName %in% commonS1S2Singlet)
finalUMAPJointSingletAll = finalUMAPJoint %>% filter(barcodeName %in% commonS1S2SingletAll)
finalUMAPJointNonSingletSinglet = finalUMAPJoint %>% filter(barcodeName %in% commonS1S2NonSingletSinglet)

#####Supplementary Figure 
plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 0.3, shape = 16) + 
  geom_point(data= finalUMAPJointSinglet, aes(UMAP_1, UMAP_2, color = sampleNum), size = 0.3, shape = 16) + 
  facet_wrap(~barcodeName, ncol = 4) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_memorySingletS1S2MLANA_individual.svg'), width = 6, height = 3.05)


plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 0.3, shape = 16) + 
  geom_point(data= finalUMAPJointSingletAll, aes(UMAP_1, UMAP_2, color = sampleNum), size = 0.3, shape = 16) + 
  facet_wrap(~barcodeName, ncol = 6) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_memorySingletS1S2MLANAALL_individual.svg'), width = 6, height = 6)

plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 0.3, shape = 16) + 
  geom_point(data= finalUMAPJointNonSingletSinglet, aes(UMAP_1, UMAP_2, color = sampleNum), size = 0.3, shape = 16) + 
  facet_wrap(~barcodeName, ncol = 4) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_memorySingletNonSingletS1S2MLANA_individual.svg'), width = 6, height = 4.21) ###thid plot likely won't need for final paper

#####density plot
plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 1, shape = 16) + 
  geom_density_2d(data= finalUMAPJointSinglet, aes(UMAP_1, UMAP_2, color = sampleNum), contour_var = "ndensity") +
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM02_memoryClusterMLANAAllSinglets.svg'), width = 6, height = 5.321)



#######Complementary analysis from S2 to S1

finalUMAPS1 = finalUMAPJointMLANA %>% filter(sampleNum=="S1")
finalUMAPS2 = finalUMAPJointMLANA %>% filter(sampleNum=="S2")
finalUMAPS2All = finalUMAPJoint %>% filter(sampleNum=="S2")
finalUMAPS1All = finalUMAPJoint %>% filter(sampleNum=="S1")

finalUMAPS2Singlet = finalUMAPS2 %>% group_by(barcodeName,sampleNum) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony ==1)
finalUMAPS1SingletSibling = finalUMAPS1 %>% group_by(barcodeName,sampleNum) %>% summarise(nColony = length(barcodeName)) 
finalUMAPS1SingletSiblingAll = finalUMAPS1All %>% group_by(barcodeName,sampleNum) %>% summarise(nColony = length(barcodeName)) 

#finalUMAPS2NonSingletSibling = finalUMAPS2 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony >1) %>% select(-nColony)


commonS1S2Singlet = inner_join(finalUMAPS2Singlet, finalUMAPS1SingletSibling, by = "barcodeName")
commonS1S2Singlet = commonS1S2Singlet$barcodeName

commonS1S2SingletAll = inner_join(finalUMAPS2Singlet, finalUMAPS1SingletSiblingAll, by = "barcodeName")
commonS1S2SingletAll = commonS1S2SingletAll$barcodeName

finalUMAPJointSinglet = finalUMAPJoint %>% filter(barcodeName %in% commonS1S2Singlet)
finalUMAPJointSingletAll = finalUMAPJoint %>% filter(barcodeName %in% commonS1S2SingletAll)

plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 0.3, shape = 16) + 
  geom_point(data= finalUMAPJointSinglet, aes(UMAP_1, UMAP_2, color = sampleNum), size = 0.3, shape = 16) + 
  facet_wrap(~barcodeName, ncol = 4) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_memorySingletS2S1MLANA_individual.svg'), width = 6, height = 3.05) ###likely won't use for paper.

plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 0.3, shape = 16) + 
  geom_point(data= finalUMAPJointSingletAll, aes(UMAP_1, UMAP_2, color = sampleNum), size = 0.3, shape = 16) + 
  facet_wrap(~barcodeName, ncol = 6) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_memorySingletS2S1MLANAALL_individual.svg'), width = 6, height = 6) ###likely won't use for paper.

#############################################################################
############################### Sibling NGFR analysis High Vs Low ######################
############20210425###########################################################
clusterNGFRCellIDs = umapClusters %>% filter(scanorama_snn_res.0.8 == 9)
clusterNGFRUMAP = inner_join(clusterNGFRCellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterNGFRUMAPBarcodes = inner_join(linCountTooverlaps, clusterNGFRUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,scanorama_snn_res.0.8)
colonySizeClusterNGFR = clusterNGFRUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterNGFR,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterNGFR = nColony.x/nColony.y) %>% filter(fractionClusterNGFR >0.5) %>% select(-nColony.x)
colonySizeClusterNGFRSummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))


#####basic plots of cluster 2 and Sample 1 vs Sample 2 
plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1, shape = 16) +
  geom_point(data = clusterNGFRUMAP, aes(x = UMAP_1, y = UMAP_2), color = "black", size = 1, shape = 16) +
  create_lpr_theme()
ggsave(plot, file = paste0(plot3Directory, 'ClusterNGFR_S1S2.svg'), width = 6, height = 5.321)

####Singlets and large colonies in S1 and S2 in cluster NGFR
colonyS1 = fractionColony %>% filter(nColony.y >1, sampleNum == "S1") 
colonyS1UMAP = jointUMAP %>% filter(sampleNum =="S1") %>% filter(BC50StarcodeD8 %in% c(colonyS1$BC50StarcodeD8))

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1, shape = 16) +
  geom_point(data = colonyS1UMAP, aes(x = UMAP_1, y = UMAP_2, col=BC50StarcodeD8), size = 2, shape = 16) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'clusterNGFRlargeColoniesS1.svg'), width = 6, height =  5.321)

colonyS2 = fractionColony %>% filter(nColony.y >1, sampleNum == "S2") 
colonyS2UMAP = jointUMAP %>% filter(sampleNum =="S2") %>% filter(BC50StarcodeD8 %in% c(colonyS2$BC50StarcodeD8))
###there is no colony >1 for S2.

SingletS1 = fractionColony %>% filter(nColony.y ==1, sampleNum == "S1") 
selectS1UMAP = jointUMAP %>% filter(sampleNum =="S1") %>% filter(BC50StarcodeD8 %in% c(SingletS1$BC50StarcodeD8))

SingletS2 = fractionColony %>% filter(nColony.y ==1, sampleNum == "S2") 
selectS2UMAP = jointUMAP %>% filter(sampleNum =="S2") %>% filter(BC50StarcodeD8 %in% c(SingletS2$BC50StarcodeD8))

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1, shape = 16) +
  geom_point(data = selectS1UMAP, aes(x = UMAP_1, y = UMAP_2, col=BC50StarcodeD8), size = 2, shape = 16) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'clusterNGFRSingletS1.svg'), width = 6, height =  5.321)

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1, shape = 16) +
  geom_point(data = selectS2UMAP, aes(x = UMAP_1, y = UMAP_2, col=BC50StarcodeD8), size = 2, shape = 16) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'clusterNGFRSingletS2.svg'), width = 6, height =  5.321)


#############################################################################
##############Sibling analysis for singlets and colonies in Sample 1 (cluster 9, NGFR), other clusters #############
#############################################################################
clusterNGFRCellIDs = umapClusters %>% filter(scanorama_snn_res.0.8 == 9)
clusterNGFRUMAP = inner_join(clusterNGFRCellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterNGFRUMAPBarcodes = inner_join(linCountTooverlaps, clusterNGFRUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,scanorama_snn_res.0.8)
colonySizeClusterNGFR = clusterNGFRUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterNGFR,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterNGFR = nColony.x/nColony.y) %>% filter(fractionClusterNGFR >0.5) %>% select(-nColony.x)
colonySizeClusterNGFRSummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

jointUMAPS1 = jointUMAP %>% filter(sampleNum=="S1")
jointUMAPS2 = jointUMAP %>% filter(sampleNum=="S2")

jointUMAPS1Barcodes = jointUMAP %>% filter(sampleNum=="S1") %>% select(BC50StarcodeD8) %>% unique()
jointUMAPS2Barcodes = jointUMAP %>% filter(sampleNum=="S2")  %>% select(BC50StarcodeD8) %>% unique()

jointBarcodesOnlyBoth = inner_join(jointUMAPS1Barcodes,jointUMAPS2Barcodes, by = "BC50StarcodeD8") 
jointBarcodesOnlyBoth = jointBarcodesOnlyBoth %>% mutate(barcodeName = c(1:nrow(jointBarcodesOnlyBoth)))
jointBarcodesOnlyBoth$barcodeName = sub("^", "B", jointBarcodesOnlyBoth$barcodeName)
jointBarcodesOnlyBothList = jointBarcodesOnlyBoth$BC50StarcodeD8

jointUMAPOnlyBoth = jointUMAP %>% filter(BC50StarcodeD8 %in% jointBarcodesOnlyBothList)
finalUMAPJoint = inner_join(jointUMAPOnlyBoth,jointBarcodesOnlyBoth, by = "BC50StarcodeD8") %>% select(-cellID, -nUMI)
finalUMAPJointClusterNGFR = inner_join(finalUMAPJoint, fractionColony, by = c("BC50StarcodeD8","sampleNum")) %>% select(-nColony.y, -fractionClusterNGFR)

finalUMAPS1 = finalUMAPJointClusterNGFR %>% filter(sampleNum=="S1")
finalUMAPS2 = finalUMAPJoint %>% filter(sampleNum=="S2")
finalUMAPS1SingletColony = finalUMAPS1 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony >0) %>% select(-nColony)
finalUMAPS2SingletNonSingletSibling = finalUMAPS2 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) 

commonS1S2SingletNonSinglet = inner_join(finalUMAPS1SingletColony, finalUMAPS2SingletNonSingletSibling)
commonS1S2SingletNonSinglet = commonS1S2SingletNonSinglet$barcodeName

finalUMAPJointSingletNonSinglet = finalUMAPJoint %>% filter(barcodeName %in% commonS1S2SingletNonSinglet)

##test similar approach to FM02_s1s3
#fractionColonyNGFRS1 = fractionColony %>% filter(sampleNum == "S1")
#barcodesNGFRS1 = fractionColonyNGFRS1$BC50StarcodeD8 %>% unique()
#umapNGFRS1AllS2 = finalUMAPJoint %>% filter(BC50StarcodeD8 %in% barcodesNGFRS1)


plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93",size = 0.3, shape = 16) + 
  geom_point(data= finalUMAPJointSingletNonSinglet, aes(UMAP_1, UMAP_2, color = sampleNum),size = 0.3, shape = 16) + 
  facet_wrap(~barcodeName, ncol = 8) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_memoryClusterNGFRS1S2ALL_individual.svg'), width = 8, height = 5.56)

#####density plot
plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93",size = 1, shape = 16) + 
  geom_density_2d(data= finalUMAPJointSingletNonSinglet, aes(UMAP_1, UMAP_2, color = sampleNum), contour_var = "ndensity") +
  #scale_color_manual(values=c("hotpink3", "turquoise3")) +
  theme_classic((base_size = 36)) + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM02_memoryClusterNGFRS1S2Contour.svg'), width = 6, height = 5.321)

####PieChart
none = nrow(fractionColony %>% filter(sampleNum == "S1")) - length(commonS1S2SingletNonSinglet)
twin = length(commonS1S2SingletNonSinglet)
dataPie = data.frame(type = c("none", "twin"),
                     value = c(none/(none+twin), twin/(none+twin)))

plot = ggplot(dataPie, aes(x="", y=value, fill=type)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  theme_classic() +
  theme(axis.text.x=element_blank(), legend.position="none") +
  geom_text(aes(label = value*100), size=5)
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s2pieNGFR.svg'), width = 6, height = 6)


#############################################################################
##############Tables for Differential Analysis ######################################
#############################################################################

barcodeNameSwitch2 = c("B118",  "B138", "B140", "B173", "B54", "B81", "B86", "B163") ##those with fateSwitch type 2
barcodeSwitch1 = finalUMAPJointSingletNonSinglet %>% filter(!barcodeName %in% barcodeNameSwitch2) %>% select(-UMAP_1, -UMAP_2, -sampleNum) %>% unique()
barcodeSwitch2 = finalUMAPJointSingletNonSinglet %>% filter(barcodeName %in% barcodeNameSwitch2) %>% select(-UMAP_1, -UMAP_2, -sampleNum) %>% unique()
cellIDsForRNASeqSwitch1 = inner_join(jointUMAP,barcodeSwitch1, by = c("BC50StarcodeD8")) %>% select(cellID, barcodeName)
cellIDsForRNASeqSwitch2 = inner_join(jointUMAP,barcodeSwitch2, by = c("BC50StarcodeD8")) %>% select(cellID, barcodeName)

normalizedCountsSwitch1 = inner_join(logNormalizedCounts,cellIDsForRNASeqSwitch1, by = "cellID") %>% select(-cellID)
normalizedCountsSwitch2 = inner_join(logNormalizedCounts,cellIDsForRNASeqSwitch2, by = "cellID") %>% select(-cellID)

tallNormalizedCountsSwitch1 = normalizedCountsSwitch1 %>% gather("geneName","count", -sampleNum, -barcodeName) %>% group_by(geneName) %>% mutate(newCount = count - min(count))
tallNormalizedCountsSwitch2 = normalizedCountsSwitch2 %>% gather("geneName","count", -sampleNum, -barcodeName) %>% group_by(geneName) %>% mutate(newCount = count - min(count))

tallNormalizedCountsNGFRControl1 = normalizedCountsSwitch1 %>% gather("geneName","count", -sampleNum, -barcodeName) %>% group_by(geneName) %>% filter(sampleNum == "S1") %>% mutate(set = "switch1") %>% select(-sampleNum)
tallNormalizedCountsNGFRControl2 = normalizedCountsSwitch2 %>% gather("geneName","count", -sampleNum, -barcodeName) %>% group_by(geneName) %>% filter(sampleNum == "S1") %>% mutate(set = "switch2") %>% select(-sampleNum)
tallNormalizedCountsNGFRControl =  bind_rows(tallNormalizedCountsNGFRControl1,tallNormalizedCountsNGFRControl2) %>% group_by(geneName) %>% mutate(newCount = count - min(count))

averagetallNormalizedCountsSwitch1 = tallNormalizedCountsSwitch1 %>% group_by(sampleNum, geneName) %>% summarise(meanCount = mean(newCount))
averagetallNormalizedCountsSwitch1Filter = averagetallNormalizedCountsSwitch1 %>% group_by(geneName) %>% filter(sum(meanCount) >1.5) 

averagetallNormalizedCountsSwitch2 = tallNormalizedCountsSwitch2 %>% group_by(sampleNum, geneName) %>% summarise(meanCount = mean(newCount))
averagetallNormalizedCountsSwitch2Filter = averagetallNormalizedCountsSwitch2 %>% group_by(geneName) %>% filter(sum(meanCount) >1.5)

averagetallNormalizedCountsNGFRControl = tallNormalizedCountsNGFRControl %>% group_by(set, geneName) %>% summarise(meanCount = mean(newCount))
averagetallNormalizedCountsNGFRControlFilter = averagetallNormalizedCountsNGFRControl %>% group_by(geneName) %>% filter(sum(meanCount) >1.5)

averagetallNormalizedCountsSwitch1FilterFinal = averagetallNormalizedCountsSwitch1Filter %>%  spread(sampleNum,meanCount) %>% mutate(foldChange = S1/S2) %>% filter(foldChange>2 | foldChange<0.5)
averagetallNormalizedCountsSwitch2FilterFinal = averagetallNormalizedCountsSwitch2Filter %>%  spread(sampleNum,meanCount) %>% mutate(foldChange = S1/S2) %>% filter(foldChange>2 | foldChange<0.5)
averagetallNormalizedCountsNGFRControlFilterFinal = averagetallNormalizedCountsNGFRControlFilter %>%  spread(set,meanCount) %>% mutate(foldChange = switch1/switch2) %>% filter(foldChange>2 | foldChange<0.5)

write.table(averagetallNormalizedCountsSwitch1FilterFinal, file=paste0(plot3Directory,'averagetallNormalizedCountsSwitch1FilterFinal.tsv'), col.names = TRUE, sep='\t')
write.table(averagetallNormalizedCountsSwitch2FilterFinal, file=paste0(plot3Directory,'averagetallNormalizedCountsSwitch2FilterFinal.tsv'), col.names = TRUE, sep='\t')
write.table(averagetallNormalizedCountsNGFRControlFilterFinal, file=paste0(plot3Directory,'averagetallNormalizedCountsNGFRControlFilterFinal.tsv'), col.names = TRUE, sep='\t')

averagetallNormalizedCountsSwitch1FilterFinal = as_tibble(read.table(file = paste0(plot3Directory, "averagetallNormalizedCountsSwitch1FilterFinal.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
averagetallNormalizedCountsSwitch2FilterFinal = as_tibble(read.table(file = paste0(plot3Directory, "averagetallNormalizedCountsSwitch2FilterFinal.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
averagetallNormalizedCountsNGFRControlFilterFinal = as_tibble(read.table(file = paste0(plot3Directory, "averagetallNormalizedCountsNGFRControlFilterFinal.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))


####Add HeatMap script here
#sort data in ascending order based on foldChange values, log2 transform data
averagetallNormalizedCountsSwitch1FilterFinal <- averagetallNormalizedCountsSwitch1FilterFinal %>% 
  arrange(foldChange) %>% 
  mutate(S1 = log2(S1)) %>% 
  mutate(S2 = log2(S2))

#pheatmap takes in a matrix as input so we need to take our numerical data as a matrix
plotTib <- as.matrix(averagetallNormalizedCountsSwitch1FilterFinal[,2:3])

#take the gene names from tibble and make it the matrix row names
rownames(plotTib) <- averagetallNormalizedCountsSwitch1FilterFinal$geneName

#transpose matrix
plotTib <- t(plotTib)

#create settings for color scaling  
breaks = seq(-2,2, by = 0.05) #hard code scale from -2 to 2
col <- colorRampPalette(brewer.pal(n=10, "RdYlBu"))(length(breaks)) #set colors
col <- rev(col)

#plot
p1 <- pheatmap(plotTib, cluster_rows = F, cluster_cols = F,
               gaps_col = 43, angle_col =45, fontsize_col = 7,
               color = col,
               breaks = breaks)

ggsave(p1, file = paste0(plot3Directory, 'FM02_pHeatMapSwitch1.svg'), width = 9, height = 3)
ggsave(p1, file = paste0(plot3Directory, 'FM02_pHeatMapSwitch1.pdf'), width = 15, height = 3)


###switch 2
averagetallNormalizedCountsSwitch2FilterFinal <- averagetallNormalizedCountsSwitch2FilterFinal %>% 
  arrange(foldChange) %>% 
  mutate(S1 = log2(S1)) %>% 
  mutate(S2 = log2(S2))

#pheatmap takes in a matrix as input so we need to take our numerical data as a matrix
plot2Tib <- as.matrix(averagetallNormalizedCountsSwitch2FilterFinal[,2:3])

#take the gene names from tibble and make it the matrix row names
rownames(plot2Tib) <- averagetallNormalizedCountsSwitch2FilterFinal$geneName

#transpose matrix
plot2Tib <- t(plot2Tib)

#create settings for color scaling  
breaks = seq(-3,3, by = 0.05) #hard code scale from -2 to 2
col <- colorRampPalette(brewer.pal(n=10, "RdYlBu"))(length(breaks)) #set colors
col <- rev(col)
#plot
p1 <- pheatmap(plot2Tib, cluster_rows = F, cluster_cols = F,
               angle_col = 45, fontsize_col = 7,
               color = col,
               breaks = breaks)

ggsave(p1, file = paste0(plot3Directory, 'FM02_pHeatMapSwitch2.svg'), width = 9, height = 3)
ggsave(p1, file = paste0(plot3Directory, 'FM02_pHeatMapSwitch2.pdf'), width = 20, height = 3)

######controlNGFR
averagetallNormalizedCountsNGFRControlFilterFinal <- averagetallNormalizedCountsNGFRControlFilterFinal %>% 
  arrange(foldChange) %>% 
  mutate(switch1 = log2(switch1)) %>% 
  mutate(switch2 = log2(switch2))

#pheatmap takes in a matrix as input so we need to take our numerical data as a matrix
plot2Tib <- as.matrix(averagetallNormalizedCountsNGFRControlFilterFinal[,2:3])

#take the gene names from tibble and make it the matrix row names
rownames(plot2Tib) <- averagetallNormalizedCountsNGFRControlFilterFinal$geneName

#transpose matrix
plot2Tib <- t(plot2Tib)

#create settings for color scaling  
breaks = seq(-2,2, by = 0.05) #hard code scale from -2 to 2
col <- colorRampPalette(brewer.pal(n=10, "RdYlBu"))(length(breaks)) #set colors
col <- rev(col)
#plot
p1 <- pheatmap(plot2Tib, cluster_rows = F, cluster_cols = F,
               angle_col = 45, fontsize_col = 7,
               color = col,
               breaks = breaks)
ggsave(p1, file = paste0(plot3Directory, 'FM02_pHeatMapNGFRControl.svg'), width = 7, height = 3)

################################################################################################################################
################################################All conditions Together##################################################
################################################################################################################################
averagetallNormalizedCountsSwitch1FilterFinal = as_tibble(read.table(file = paste0(plot3Directory, "averagetallNormalizedCountsSwitch1FilterFinal.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
averagetallNormalizedCountsSwitch2FilterFinal = as_tibble(read.table(file = paste0(plot3Directory, "averagetallNormalizedCountsSwitch2FilterFinal.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
averagetallNormalizedCountsNGFRControlFilterFinal = as_tibble(read.table(file = paste0(plot3Directory, "averagetallNormalizedCountsNGFRControlFilterFinal.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))

geneListAll = bind_rows(averagetallNormalizedCountsSwitch1FilterFinal, averagetallNormalizedCountsSwitch2FilterFinal, averagetallNormalizedCountsNGFRControlFilterFinal) %>% select(geneName) %>% unique()

averagetallNormalizedCountsSwitch1Compare = averagetallNormalizedCountsSwitch1 %>% filter (geneName %in% geneListAll$geneName)
averagetallNormalizedCountsSwitch2Compare = averagetallNormalizedCountsSwitch2 %>% filter (geneName %in% geneListAll$geneName)

ngfrSwitch1 = averagetallNormalizedCountsSwitch1Compare %>% filter(sampleNum == "S1") %>% mutate(type = "ngfrSwitch1")
ngfrSwitch2 = averagetallNormalizedCountsSwitch2Compare %>% filter(sampleNum == "S1") %>% mutate(type = "ngfrSwitch2")
switch1 = averagetallNormalizedCountsSwitch1Compare %>% filter(sampleNum == "S2") %>% mutate(type = "switch1")
switch2 = averagetallNormalizedCountsSwitch2Compare %>% filter(sampleNum == "S2") %>% mutate(type = "switch2")

allTable = bind_rows(ngfrSwitch1, ngfrSwitch2, switch1, switch2) %>% ungroup() %>% select(-sampleNum)
allTableSpreadTibble = allTable %>% spread(type, meanCount)
allTableSpread = as.matrix(allTableSpreadTibble[,2:5])
rownames(allTableSpread) = allTableSpreadTibble$geneName
allTableSpread <- t(allTableSpread)

#create settings for color scaling  
breaks = seq(-2,2, by = 0.05) #hard code scale from -2 to 2
col <- colorRampPalette(brewer.pal(n=10, "RdYlBu"))(length(breaks)) #set colors
col <- rev(col)
#plot
p1 <- pheatmap(allTableSpread, cluster_rows = T, cluster_cols = T,
               angle_col = 90, fontsize_col = 4,fontsize_row = 4, scale = 'column',
               color = col)
ggsave(p1, file = paste0(plot3Directory, 'FM02_AllComparison.svg'), width = 14, height = 4)

################################################################################################################################


################################################################################################################################
################################################Baseline Differential analysis##################################################
################################################################################################################################
###select random 150 from each
finalUMAPJointSingletNonSingletS1 = finalUMAPJointSingletNonSinglet %>% filter(sampleNum == "S1")
logNormalizedCountsS1 = logNormalizedCounts %>% filter(sampleNum == "S1")
randomSampling = c()
randomSampling2 = c()

for (i in 1:1000) {
  
  n1=sample(15:100,1)
  n2=sample(15:100,1)
  
  a = tibble(barcodeSampling1$UMAP_1) %>% dplyr::rename(value = "barcodeSampling1$UMAP_1")
  b = tibble(barcodeSampling2$UMAP_1) %>% dplyr::rename(value = "barcodeSampling2$UMAP_1")
  
  barcodeSampling1 = finalUMAPJointSingletNonSingletS1[sample(nrow(finalUMAPJointSingletNonSingletS1), n1), ] %>% select(-sampleNum)
  finalUMAPJointSingletNonSingletS1Rep = finalUMAPJointSingletNonSingletS1 %>% filter(!UMAP_1 %in% c(barcodeSampling1$UMAP_1))
  barcodeSampling2 = finalUMAPJointSingletNonSingletS1Rep[sample(nrow(finalUMAPJointSingletNonSingletS1Rep), n2), ] %>% select(-sampleNum)
  
  cellIDsForRNASeqSampling1 = inner_join(jointUMAP,barcodeSampling1, by = c("BC50StarcodeD8", "UMAP_1", "UMAP_2")) %>% select(cellID, barcodeName)
  cellIDsForRNASeqSampling2 = inner_join(jointUMAP,barcodeSampling2, by = c("BC50StarcodeD8","UMAP_1", "UMAP_2")) %>% select(cellID, barcodeName)
  
  normalizedCountsSampling1 = inner_join(logNormalizedCountsS1,cellIDsForRNASeqSampling1, by = "cellID") %>% select(-cellID)
  normalizedCountsSampling2 = inner_join(logNormalizedCountsS1,cellIDsForRNASeqSampling2, by = "cellID") %>% select(-cellID)
  
  tallNormalizedCountsSampling1 = normalizedCountsSampling1 %>% gather("geneName","count", -sampleNum, -barcodeName) %>% group_by(geneName) %>% mutate(set = "switch1") %>% select(-sampleNum)
  tallNormalizedCountsSampling2 = normalizedCountsSampling2 %>% gather("geneName","count", -sampleNum, -barcodeName) %>% group_by(geneName) %>% mutate(set = "switch2") %>% select(-sampleNum)
  tallNormalizedCountsNGFRSampling =  bind_rows(tallNormalizedCountsSampling1,tallNormalizedCountsSampling2) %>% group_by(geneName) %>% mutate(newCount = count - min(count))
  
  averagetallNormalizedCountsNGFRSampling = tallNormalizedCountsNGFRSampling %>% group_by(set, geneName) %>% summarise(meanCount = mean(newCount))
  averagetallNormalizedCountsNGFRSamplingFilter = averagetallNormalizedCountsNGFRSampling %>% group_by(geneName) %>% filter(sum(meanCount) >1.5)
  
  averagetallNormalizedCountsNGFRSamplingFilterFinal = averagetallNormalizedCountsNGFRSamplingFilter %>%  spread(set,meanCount) %>% mutate(foldChange = switch1/switch2) %>% filter(foldChange>2 | foldChange<0.5)
  
  randomSampling2[i] = nrow(averagetallNormalizedCountsNGFRSamplingFilterFinal)
}

randomSampling1 = as_tibble(randomSampling)
randomSampling2 = as_tibble(randomSampling2)

write.table(randomSampling2, file=paste0(plot3Directory,'FM02_randomSamplingNGFR2.tsv'), col.names = TRUE, sep='\t')

mean(randomSampling2$value) # 20.44
median(randomSampling2$value) # 11

pvalueNGFR = length(randomSampling2[randomSampling2>(nrow(averagetallNormalizedCountsNGFRControlFilterFinal)-1)])/nrow(randomSampling2)
pvalueswitch1 = length(randomSampling2[randomSampling2>(nrow(averagetallNormalizedCountsSwitch1FilterFinal)-1)])/nrow(randomSampling2)
pvalueswitch2 = length(randomSampling2[randomSampling2>(nrow(averagetallNormalizedCountsSwitch2FilterFinal)-1)])/nrow(randomSampling2)



plot = ggplot(randomSampling1, aes(x=value)) + geom_histogram() + 
  geom_vline(aes(xintercept=mean(value)),color="blue", size=1) +
  geom_vline(aes(xintercept=nrow(averagetallNormalizedCountsNGFRControlFilterFinal)),color="blue", size=1) +
  geom_vline(aes(xintercept=nrow(averagetallNormalizedCountsSwitch2FilterFinal)),color="blue", size=1) +
  geom_vline(aes(xintercept=nrow(averagetallNormalizedCountsSwitch1FilterFinal)),color="blue", size=1) +
  theme_classic()
  
ggsave(plot, file = paste0(plot3Directory, 'FM02_randomSamplingWithReplacement_NGFR.svg'), width = 4, height = 4)

plot = ggplot(randomSampling2, aes(x=value)) + geom_histogram(binwidth=5) +theme_classic() + 
  geom_vline(aes(xintercept=mean(value)),color="blue", size=1) +
  geom_vline(aes(xintercept=nrow(averagetallNormalizedCountsNGFRControlFilterFinal)),color="blue", size=1) +
  geom_vline(aes(xintercept=nrow(averagetallNormalizedCountsSwitch2FilterFinal)),color="blue", size=1) +
  geom_vline(aes(xintercept=nrow(averagetallNormalizedCountsSwitch1FilterFinal)),color="blue", size=1)+
  theme_classic()
ggsave(plot, file = paste0(plot3Directory, 'FM02_randomSamplingWithoutReplacement_NGFR.svg'), width = 5, height = 4)

########################################################################################################################################
###########################################BARCODES FOR 10X BARCODE COMPARISONS##############################################
########################################################################################################################################

lowDoseDependentBarcodesReverseComplement = as_tibble(read.table(paste0(home3Directory, 'lowDoseDependentBarcodesReverseComplement.csv'), stringsAsFactors=F, header = F))
lowDoseDependentBarcodesReverseComplement = lowDoseDependentBarcodesReverseComplement %>% mutate(barcodeSub = substr(lowDoseDependentBarcodesReverseComplement$V1,1,50)) %>% select(-V1)
nrow(stringdist_inner_join(lowDoseDependentBarcodesReverseComplement,linCountTooverlaps,by=c('barcodeSub'='BC50StarcodeD8'),max_dist=4)) ###this one works

overlapBarcodes = stringdist_inner_join(lowDoseDependentBarcodesReverseComplement,linCountTooverlaps,by=c('barcodeSub'='BC50StarcodeD8'),max_dist=4)
overlapBarcodesS2 = overlapBarcodes %>% filter(sampleNum == "S2") %>% select (cellID, BC50StarcodeD8) %>% unique()

S2AllCoordinates = jointUMAP %>% filter(sampleNum == "S2")
S2OnlyCoordinates = inner_join(overlapBarcodesS2, umapCoordinates, by = c("cellID"))

plot <- ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93",size = 1, shape = 16) +
  geom_point(data = S2AllCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3",size = 2, shape = 16) +
  geom_point(data = S2OnlyCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "green3",size = 2, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plot3Directory, 'FM02_AllColoniesOverlayLowDoseAndLowDoseOnly.svg'), width = 6, height = 5.321)


########################## PLOT for gDNASeq analysis #####################
ratio_highLow = c(0.6193759943,0.4296939031, 0.5798988486)
ratio_highHigh = c(1.029762575,1.073929232,0.9522726901)

highLow = tibble(ratio = numeric()) %>% add_row(ratio = ratio_highLow) %>% mutate(type = "highLow")
highHigh = tibble(ratio = numeric()) %>% add_row(ratio = ratio_highHigh) %>% mutate(type = "highHigh")

all = bind_rows(highLow, highHigh)
plot2 = ggplot(all, aes(x = type, y=ratio, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  geom_jitter() +
  scale_fill_manual(values=c("grey", "grey")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none") + 
  ylim(0,1.2)

ggsave(plot2, file = paste0(plot5Directory, 'UniqueBarcodes.svg'), width = 4, height = 8) ###not sure whether we will use it.

