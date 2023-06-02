#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2020_FM03/jointAnalysis/seurat/hg19/DOT1Lsample1Sample2_V1_hg19.R
##############################################

library(Seurat)
library(ggplot2)
library(tidyverse)
library(fuzzyjoin)
options(future.globals.maxSize = 4000 * 1024^2)

#########
# workflow adapted from https://satijalab.org/seurat/v3.0/integration.html
#########
home1Directory <-"/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM03/"
home2Directory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM03/'
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

s1s2Scanorama <- readRDS(paste0(home1Directory,"s1s2_fm03_scanorama_50pcs_filter_MT30_hg19.rds"))
logNormalizedCounts = as_tibble(read.table(file = paste0(home1Directory, "logNormalizedCounts_s1s2Scanorama_50pcs_filterRound_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinates = as_tibble(read.table(file = paste0(home1Directory, "umapCoordinates_s1s2Scanorama_50pcs_filter_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapClusters = as_tibble(read.table(file = paste0(home1Directory, "umapClusters_s1s2Scanorama_50pcs_filter_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))

################Figure################################################
#####Just to plot the three with the right dimensions
s1s2Scanorama <- FindClusters(object=s1s2Scanorama,graph.name = "scanorama_snn", resolution=0.8)

clusterUMAP_08 = DimPlot(s1s2Scanorama, reduction = "umap_scanorama", label = FALSE, group.by = "seurat_clusters") + NoLegend()
clusterUMAP_08WCluster = DimPlot(s1s2Scanorama, reduction = "umap_scanorama", label = TRUE, label.size = 8, group.by = "seurat_clusters") + NoLegend()
sampleUMAP = DimPlot(s1s2Scanorama, reduction = "umap_scanorama", group.by = "orig.ident") + NoLegend()
sampleSplitUMAP_08 = DimPlot(s1s2Scanorama, reduction = "umap_scanorama",slot="RNA.data", group.by = "seurat_clusters", split.by = "orig.ident", label = FALSE) + NoLegend()

ggsave(clusterUMAP_08, file = paste0(plot3Directory, 'FM03_s1s2_clusterUMAP_onlyScanorma08_hg19.svg'), width = 6, height = 5.321)
ggsave(clusterUMAP_08WCluster, file = paste0(plot3Directory, 'FM03_s1s2_clusterUMAPWNUMBERS_onlyScanorma08_hg19.svg'), width = 6, height = 5.321)
ggsave(sampleUMAP, file = paste0(plot3Directory, 'FM03_s1s2_sampleUMAP_onlyScanorma_hg19.svg'), width = 6, height = 5.321)
ggsave(sampleSplitUMAP_08, file = paste0(plot3Directory, 'FM03_s1s2_sampleSplitUMAP_onlyScanorma08_hg19.svg'), width = 12, height = 5.822)

plot <- ggplot(umapCoordinates, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes_string(color = logNormalizedCounts$APOE), size = 1, shape = 16) +
  scale_color_gradient(low = "gray93", high = "darkblue") +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Scaled\n Log UMI counts")
ggsave(plot, file = paste0(plot3Directory, 'FM03_MLANA.svg'), width = 6, height = 6.429)

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

jointUMAPS1 = jointUMAP %>% filter(sampleNum == "S1")
jointUMAPS2 = jointUMAP %>% filter(sampleNum == "S2")

length(unique(jointUMAPS2$BC50StarcodeD8))
#############################################################################
##############Analysis for MLANA Cluster High Vs Low ########################
#############################################################################
#nDMSO = 9590
#nDOT1L = 7547 
#nRatio = nDMSO/nDOT1L

clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(2) )
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))

colonySizeAllToPlot = colonySizeAll %>% mutate(type = case_when(nColony == 1 ~ "singlet", 
                                                                nColony == 2 ~ "doublet",
                                                                nColony > 2 ~ "colony"),
                                               type2 = case_when(nColony == 1 ~ "singlet", 
                                                                 nColony >1 ~ "colony"))
colonySizeAllToPlot2 = colonySizeAllToPlot %>% group_by(sampleNum, type) %>% summarise(aggreate = length(BC50StarcodeD8)) %>% 
  ungroup() %>% group_by(sampleNum) %>%mutate(fraction = aggreate/sum(aggreate))
colonySizeAllToPlot2$type <- factor(colonySizeAllToPlot2$type, levels=c("colony", "doublet", "singlet"))
colonySizeAllToPlot2$sampleNum <- factor(colonySizeAllToPlot2$sampleNum, levels=c("S2", "S1"))

p <- ggplot(data=colonySizeAllToPlot2, aes(x=type, y=fraction, fill=sampleNum)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal() + coord_flip() +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())

colonySizeAllToPlot3 = colonySizeAllToPlot %>% group_by(sampleNum, type2) %>% summarise(aggreate = length(BC50StarcodeD8)) %>% 
  ungroup() %>% group_by(sampleNum) %>%mutate(fraction = aggreate/sum(aggreate))
colonySizeAllToPlot3$type2 <- factor(colonySizeAllToPlot3$type2, levels=c("colony", "singlet"))
colonySizeAllToPlot3$sampleNum <- factor(colonySizeAllToPlot3$sampleNum, levels=c("S2", "S1"))

plot <- ggplot(data=colonySizeAllToPlot3, aes(x=type2, y=fraction, fill=sampleNum)) +
  geom_bar(stat="identity", color="black", position=position_dodge(), width=0.3)+
  theme_minimal() + coord_flip() +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM03_singletColonySummary.svg'), width = 6, height = 4)

######################################


fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

#####basic plots of cluster 2 and Sample 1 vs Sample 2 

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1, shape = 16) +
  geom_point(data = clusterMLANAUMAP, aes(x = UMAP_1, y = UMAP_2, color = factor(seurat_clusters)), size = 1, shape = 16) +
  create_lpr_theme()+ theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM03_ClusterMLANA_S1S2.svg'), width = 6, height = 5.321)

colonyS2 = fractionColony %>% filter(nColony.y >2, sampleNum == "S2") 
colonyS2UMAP = jointUMAP %>% filter(sampleNum =="S2") %>% filter(BC50StarcodeD8 %in% c(colonyS2$BC50StarcodeD8))

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93",size = 2, shape = 16) +
  geom_point(data = colonyS2UMAP, aes(x = UMAP_1, y = UMAP_2, col=BC50StarcodeD8),size = 3, shape = 16) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM03_largeColoniesS2MLANA.svg'), width = 6, height = 5.321)

SingletS1 = fractionColony %>% filter(nColony.y ==1, sampleNum == "S1") 
selectS1UMAP = jointUMAP %>% filter(sampleNum =="S1") %>% filter(BC50StarcodeD8 %in% c(SingletS1$BC50StarcodeD8))

SingletS1All = colonySizeAll %>% filter(nColony ==1, sampleNum == "S1") ####just to check the distributions
selectS1UMAPAll = jointUMAP %>% filter(sampleNum =="S1") %>% filter(BC50StarcodeD8 %in% c(SingletS1All$BC50StarcodeD8)) ####just to check the distributions

SingletS2All = colonySizeAll %>% filter(nColony ==1, sampleNum == "S2") ####just to check the distributions
selectS2UMAPAll = jointUMAP %>% filter(sampleNum =="S2") %>% filter(BC50StarcodeD8 %in% c(SingletS2All$BC50StarcodeD8)) ####just to check the distributions

SingletS2 = fractionColony %>% filter(nColony.y ==1, sampleNum == "S2") 
selectS2UMAP = jointUMAP %>% filter(sampleNum =="S2") %>% filter(BC50StarcodeD8 %in% c(SingletS2$BC50StarcodeD8))


plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93",size = 2, shape = 16) +
  geom_point(data = selectS1UMAP, aes(x = UMAP_1, y = UMAP_2), col="hotpink3", size = 3, shape = 16) +
  geom_point(data = selectS2UMAP, aes(x = UMAP_1, y = UMAP_2), col="turquoise3",size = 3, shape = 16) +
  create_lpr_theme()
ggsave(plot, file = paste0(plot3Directory, 'FM03_SingletsS1S2_ClusterMLANA_C2.svg'), width = 6, height = 5.321)


plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93",size = 2, shape = 16) +
  geom_point(data = selectS1UMAPAll, aes(x = UMAP_1, y = UMAP_2), col="hotpink3", size = 3, shape = 16) +
  geom_point(data = selectS2UMAPAll, aes(x = UMAP_1, y = UMAP_2), col="turquoise3",size = 3, shape = 16) +
  create_lpr_theme()
ggsave(plot, file = paste0(plot3Directory, 'FM03_SingletsS1S2_All.svg'), width = 6, height = 5.321)

###added on 10/10/2021 for colony
ColonyS1All = colonySizeAll %>% filter(nColony >1, sampleNum == "S1") ####just to check the distributions
colonyS1UMAPAll = jointUMAP %>% filter(sampleNum =="S1") %>% filter(BC50StarcodeD8 %in% c(ColonyS1All$BC50StarcodeD8)) ####just to check the distributions

ColonyS2All = colonySizeAll %>% filter(nColony >1, sampleNum == "S2") ####just to check the distributions
colonyS2UMAPAll = jointUMAP %>% filter(sampleNum =="S2") %>% filter(BC50StarcodeD8 %in% c(ColonyS2All$BC50StarcodeD8)) ####just to check the distributions

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93",size = 2, shape = 16) +
  geom_point(data = colonyS1UMAPAll, aes(x = UMAP_1, y = UMAP_2), col="hotpink3", size = 3, shape = 16) +
  geom_point(data = colonyS2UMAPAll, aes(x = UMAP_1, y = UMAP_2), col="turquoise3",size = 3, shape = 16) +
  create_lpr_theme()
ggsave(plot, file = paste0(plot3Directory, 'FM03_ColonyS1S2_All.svg'), width = 6, height = 5.321)

#####

#############################################################################
####### singlet ONLY Cluster MLANA ######################
#############################################################################
clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(2) )
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
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
finalUMAPJointMLANA = inner_join(finalUMAPJoint, fractionColony, by = c("BC50StarcodeD8","sampleNum")) %>% select(-nColony.y)

finalUMAPS1 = finalUMAPJointMLANA %>% filter(sampleNum=="S1")
finalUMAPS2 = finalUMAPJointMLANA %>% filter(sampleNum=="S2")
finalUMAPS2All = finalUMAPJoint %>% filter(sampleNum=="S2")
finalUMAPS1SingletOthers = finalUMAPS1 %>% group_by(barcodeName,sampleNum) %>% summarise(nColony = length(barcodeName))
finalUMAPS2SingletSibling = finalUMAPS2 %>% group_by(barcodeName,sampleNum) %>% summarise(nColony = length(barcodeName)) 
finalUMAPS2SingletSiblingAll = finalUMAPS2All %>% group_by(barcodeName,sampleNum) %>% summarise(nColony = length(barcodeName)) 

commonS1S2SingletAll = inner_join(finalUMAPS1SingletOthers, finalUMAPS2SingletSibling, by = "barcodeName")
commonS1S2SingletAll = commonS1S2SingletAll$barcodeName

commonS1S2SingletNonSinglet = inner_join(finalUMAPS1SingletOthers, finalUMAPS2SingletSiblingAll, by = "barcodeName")
commonS1S2NonSingletSinglet = commonS1S2SingletNonSinglet$barcodeName

finalUMAPJointSingletAll = finalUMAPJoint %>% filter(barcodeName %in% commonS1S2SingletAll)
finalUMAPJointNonSingletSinglet = finalUMAPJoint %>% filter(barcodeName %in% commonS1S2NonSingletSinglet)

#####Memory Figures for DOT1L
plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 0.3, shape = 16) + 
  geom_point(data= finalUMAPJointSingletAll, aes(UMAP_1, UMAP_2, color = sampleNum), size = 0.3, shape = 16) + 
  facet_wrap(~barcodeName, ncol = 6) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM03_memorySingletS1S2MLANAALL_individual.svg'), width = 6, height = 5.321)

plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 0.3, shape = 16) + 
  geom_point(data= finalUMAPJointNonSingletSinglet, aes(UMAP_1, UMAP_2, color = sampleNum), size = 0.3, shape = 16) + 
  facet_wrap(~barcodeName, ncol = 4) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM03_memorySingletNonSingletS1S2MLANA_individual.svg'), width = 6, height = 5.8)

#####density plot
plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 1, shape = 16) + 
  geom_density_2d(data= finalUMAPJointNonSingletSinglet, aes(UMAP_1, UMAP_2, color = sampleNum), contour_var = "ndensity") +
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM03_memoryClusterMLANAAllSinglets.svg'), width = 6, height = 5.321)

####PieChart
none = nrow(fractionColony %>% filter(sampleNum == "S1")) - length(unique(finalUMAPJointNonSingletSinglet$BC50StarcodeD8))
twin = length(unique(finalUMAPJointNonSingletSinglet$BC50StarcodeD8))
dataPie = data.frame(type = c("none", "twin"),
                     value = c(none/(none+twin), twin/(none+twin)))

plot = ggplot(dataPie, aes(x="", y=value, fill=type)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  theme_classic() +
  theme(axis.text.x=element_blank(), legend.position="none") +
  geom_text(aes(label = value*100), size=5)
ggsave(plot, file = paste0(plot3Directory, 'FM03_s1s2pieMLANA.svg'), width = 6, height = 6)

#############################################################################
##############Tables for Differential Analysis ######################################
#############################################################################
barcodeSwitch = finalUMAPJointNonSingletSinglet %>% select(-UMAP_1, -UMAP_2, -sampleNum) %>% unique()
cellIDsForRNASeqSwitch = inner_join(jointUMAP,barcodeSwitch, by = c("BC50StarcodeD8")) %>% select(cellID, barcodeName)

normalizedCountsSwitch = inner_join(logNormalizedCounts,cellIDsForRNASeqSwitch, by = "cellID") %>% select(-cellID)
tallNormalizedCountsSwitch = normalizedCountsSwitch %>% gather("geneName","count", -sampleNum, -barcodeName) %>% group_by(geneName) %>% mutate(newCount = count - min(count))

averagetallNormalizedCountsSwitch = tallNormalizedCountsSwitch %>% group_by(sampleNum, geneName) %>% summarise(meanCount = mean(newCount))
averagetallNormalizedCountsSwitchFilter = averagetallNormalizedCountsSwitch %>% group_by(geneName) %>% filter(sum(meanCount) >1.5) 

averagetallNormalizedCountsSwitchFilterFinal = averagetallNormalizedCountsSwitchFilter %>%  spread(sampleNum,meanCount) %>% mutate(foldChange = S1/S2) %>% filter(foldChange>2 | foldChange<0.5)

####Add HeatMap script here
#sort data in ascending order based on foldChange values, log2 transform data
averagetallNormalizedCountsSwitchFilterFinal <- averagetallNormalizedCountsSwitchFilterFinal %>% 
  arrange(foldChange) %>% 
  mutate(S1 = log2(S1)) %>% 
  mutate(S2 = log2(S2))

#pheatmap takes in a matrix as input so we need to take our numerical data as a matrix
plotTib <- as.matrix(averagetallNormalizedCountsSwitchFilterFinal[,2:3])

#take the gene names from tibble and make it the matrix row names
rownames(plotTib) <- averagetallNormalizedCountsSwitchFilterFinal$geneName

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

ggsave(p1, file = paste0(plot3Directory, 'FM03_pHeatMapSwitch.svg'), width = 9, height = 3)
ggsave(p1, file = paste0(plot3Directory, 'FM03_pHeatMapSwitch.pdf'), width = 22, height = 3)


###########################Bar Plot for Colonies between Control and DOT1L

colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >2)

colonySizeAllS1 = colonySizeAll %>% filter(sampleNum == "S1")
colonySizeAllS2 = colonySizeAll %>% filter(sampleNum == "S2")

clustersBarcode = inner_join(umapClusters, jointUMAP, by = c("cellID","sampleNum")) %>% select(-nUMI, -UMAP_1, -UMAP_2)
clustersOnlyColoniesS1 = inner_join(clustersBarcode,colonySizeAllS1, by = c("BC50StarcodeD8","sampleNum"))
clustersOnlyColoniesS2 = inner_join(clustersBarcode,colonySizeAllS2, by = c("BC50StarcodeD8","sampleNum"))

clustersOnlyColonies = bind_rows(clustersOnlyColoniesS1,clustersOnlyColoniesS2)



plot <- ggplot(clustersOnlyColonies, aes(seurat_clusters)) +
  geom_bar(aes(y = (..count..)/sum(..count..))) +
facet_wrap(~sampleNum, ncol = 1) +
create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM03_colonyClusterCOmparisonDOT1LControl.svg'), width = 6, height = 5.321)  