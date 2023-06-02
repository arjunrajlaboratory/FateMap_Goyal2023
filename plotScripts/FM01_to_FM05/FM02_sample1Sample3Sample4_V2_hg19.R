#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2020_FM02/jointAnalysis/seurat/hg19/sample1Sample3Sample4_V1_hg19.R
##############################################

library(Seurat)
library(ggplot2)
library(tidyverse)
library(fuzzyjoin)
options(future.globals.maxSize = 4000 * 1024^2)

#########
# workflow adapted from https://satijalab.org/seurat/v3.0/integration.html
#########
home1Directory <-"/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM02/"
home2Directory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM02/'
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

s1s3s4Scanorama <- readRDS(paste0(home1Directory,"s1s3s4_scanorama_50pcs_filter_hg19.rds"))
logNormalizedCounts = as_tibble(read.table(file = paste0(home1Directory, "logNormalizedCounts_s1s3s4Scanorama_50pcs_filterRound_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinates = as_tibble(read.table(file = paste0(home1Directory, "umapCoordinates_s1s3s4Scanorama_50pcs_filter_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapClusters = as_tibble(read.table(file = paste0(home1Directory, "umapClusters_s1s3s4Scanorama_snn05_50pcs_filter_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))

################FigureBasicPlots################################################
#####Just to plot the three with the right dimensions
s1s3s4Scanorama <- FindClusters(object=s1s3s4Scanorama,graph.name = "scanorama_snn", resolution=0.5)

DimPlot(s1s3s4Scanorama, reduction = "umap_scanorama", group.by = "orig.ident") + NoLegend()
clusterUMAP_05 = DimPlot(s1s3s4Scanorama, reduction = "umap_scanorama", label = FALSE, group.by = "scanorama_snn_res.0.5") + NoLegend()
clusterUMAP_08 = DimPlot(s1s3s4Scanorama, reduction = "umap_scanorama", label = FALSE, group.by = "scanorama_snn_res.0.8") + NoLegend()
clusterUMAP_08WCluster = DimPlot(s1s3s4Scanorama, reduction = "umap_scanorama", label = TRUE,label.size = 8, group.by = "scanorama_snn_res.0.8") + NoLegend()
clusterUMAP_05WCluster = DimPlot(s1s3s4Scanorama, reduction = "umap_scanorama", label = TRUE,label.size = 8, group.by = "scanorama_snn_res.0.5") + NoLegend()
sampleUMAP = DimPlot(s1s3s4Scanorama, reduction = "umap_scanorama", group.by = "orig.ident") + NoLegend()

sampleSplitUMAP_08 = DimPlot(s1s3s4Scanorama, reduction = "umap_scanorama",slot="RNA.data", group.by = "scanorama_snn_res.0.8", split.by = "orig.ident", label = TRUE, label.size = 8,) + NoLegend()
sampleSplitUMAP_05 = DimPlot(s1s3s4Scanorama, reduction = "umap_scanorama",slot="RNA.data", group.by = "scanorama_snn_res.0.5", split.by = "orig.ident", label = TRUE, label.size = 8,) + NoLegend()

ggsave(clusterUMAP_08, file = paste0(plot3Directory, 's1s3s4_clusterUMAP_onlyScanorma08_hg19.svg'), width = 6, height = 5.321)
ggsave(clusterUMAP_05, file = paste0(plot3Directory, 's1s3s4_clusterUMAP_onlyScanorma05_hg19.svg'), width = 6, height = 5.321)
ggsave(clusterUMAP_08WCluster, file = paste0(plot3Directory, 's1s3s4_clusterUMAPWNUMBERS_onlyScanorma08_hg19.svg'), width = 6, height = 5.321)
ggsave(clusterUMAP_05WCluster, file = paste0(plot3Directory, 's1s3s4_clusterUMAPWNUMBERS_onlyScanorma05_hg19.svg'), width = 6, height = 5.321)
ggsave(sampleUMAP, file = paste0(plot3Directory, 's1s3s4_sampleUMAP_onlyScanorma_hg19.svg'), width = 6, height = 5.321)
ggsave(sampleSplitUMAP_08, file = paste0(plot3Directory, 's1s3s4_sampleSplitUMAP_onlyScanorma08WNUMBERS_hg19.svg'), width = 18, height = 5.822)
ggsave(sampleSplitUMAP_05, file = paste0(plot3Directory, 's1s3s4_sampleSplitUMAP_onlyScanorma05WNUMBERS_hg19.svg'), width = 18, height = 5.822)

plot <- ggplot(umapCoordinates, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes_string(color = logNormalizedCounts$VCAM1), size = 1, shape = 16) +
  scale_color_gradient(low = "gray93", high = "darkblue") +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Scaled\n Log UMI counts")
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3s4_VCAM1.svg'), width = 6, height = 6.429)

######implemented until here


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
  select(-c(V3,V5,V6,V7)) %>% unique() %>% select(-UMI) %>% filter(sampleNum %in% c("1","3"))

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

###ClustersToAnalyze = 3 (MLANA) and 4 (NGFR). Maybe 9?

#############################################################################
##############Analysis for MLANA Cluster High Vs Low ########################
#############################################################################

clusterMLANACellIDs = umapClusters %>% filter(scanorama_snn_res.0.5 %in% 3)
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,scanorama_snn_res.0.5)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

SingletS1All = colonySizeAll %>% filter(nColony ==1, sampleNum == "S1") ####just to check the distributions
selectS1UMAPAll = jointUMAP %>% filter(sampleNum =="S1") %>% filter(BC50StarcodeD8 %in% c(SingletS1All$BC50StarcodeD8)) ####just to check the distributions

SingletS3All = colonySizeAll %>% filter(nColony ==1, sampleNum == "S3") ####just to check the distributions
selectS3UMAPAll = jointUMAP %>% filter(sampleNum =="S3") %>% filter(BC50StarcodeD8 %in% c(SingletS1All$BC50StarcodeD8)) ####just to check the distributions

ColonyS1All = colonySizeAll %>% filter(nColony >1, sampleNum == "S1") ####just to check the distributions
ColonyS1UMAPAll = jointUMAP %>% filter(sampleNum =="S1") %>% filter(BC50StarcodeD8 %in% c(ColonyS1All$BC50StarcodeD8)) ####just to check the distributions

ColonyS3All = colonySizeAll %>% filter(nColony >1, sampleNum == "S3") ####just to check the distributions
ColonyS3UMAPAll = jointUMAP %>% filter(sampleNum =="S3") %>% filter(BC50StarcodeD8 %in% c(ColonyS3All$BC50StarcodeD8)) ####just to check the distributions

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93",size = 1, shape = 16) +
  geom_point(data = selectS1UMAPAll, aes(x = UMAP_1, y = UMAP_2), col="hotpink3", size = 2, shape = 16) +
  geom_point(data = selectS3UMAPAll, aes(x = UMAP_1, y = UMAP_2), col="turquoise3",size = 2, shape = 16) +
  create_lpr_theme()
ggsave(plot, file = paste0(plot3Directory, 'FM02_SingletsS1S3S4_All.svg'), width = 6, height = 5.321)

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93",size = 1, shape = 16) +
  geom_point(data = ColonyS1UMAPAll, aes(x = UMAP_1, y = UMAP_2), col="hotpink3", size = 2, shape = 16) +
  geom_point(data = ColonyS3UMAPAll, aes(x = UMAP_1, y = UMAP_2), col="turquoise3",size = 2, shape = 16) +
  create_lpr_theme()
ggsave(plot, file = paste0(plot3Directory, 'FM02_ColonyS1S3S4_All.svg'), width = 6, height = 5.321)

#############################################################################
#####Where do MLANA in Vemurafenib go in Trametinib
#############################################################################

clusterMLANACellIDs = umapClusters %>% filter(scanorama_snn_res.0.5 %in% 3)
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,scanorama_snn_res.0.5)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

jointUMAPS1 = jointUMAP %>% filter(sampleNum=="S1")
jointUMAPS3 = jointUMAP %>% filter(sampleNum=="S3")

jointUMAPS1Barcodes = jointUMAP %>% filter(sampleNum=="S1") %>% select(BC50StarcodeD8) %>% unique()
jointUMAPS3Barcodes = jointUMAP %>% filter(sampleNum=="S3")  %>% select(BC50StarcodeD8) %>% unique()

jointBarcodesOnlyBoth = inner_join(jointUMAPS1Barcodes,jointUMAPS3Barcodes, by = "BC50StarcodeD8") 
jointBarcodesOnlyBoth = jointBarcodesOnlyBoth %>% mutate(barcodeName = c(1:nrow(jointBarcodesOnlyBoth)))
jointBarcodesOnlyBoth$barcodeName = sub("^", "B", jointBarcodesOnlyBoth$barcodeName)
jointBarcodesOnlyBothList = jointBarcodesOnlyBoth$BC50StarcodeD8

jointUMAPOnlyBoth = jointUMAP %>% filter(BC50StarcodeD8 %in% jointBarcodesOnlyBothList)
finalUMAPJoint = inner_join(jointUMAPOnlyBoth,jointBarcodesOnlyBoth, by = "BC50StarcodeD8") %>% select(-cellID, -nUMI)
fractionColonyMLANAS1 = fractionColony %>% filter(sampleNum == "S1")
barcodesMLANAS1 = fractionColonyMLANAS1$BC50StarcodeD8 %>% unique()

umapMLANAS1AllS2 = finalUMAPJoint %>% filter(BC50StarcodeD8 %in% barcodesMLANAS1)

plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 0.3, shape = 16) + 
  geom_point(data= umapMLANAS1AllS2, aes(UMAP_1, UMAP_2, color = sampleNum), size = 0.3, shape = 16) + 
  facet_wrap(~barcodeName, ncol = 7) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3memoryS1S2MLANAALL_individual.svg'), width = 4, height = 4)

plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 1, shape = 16) + 
  geom_density_2d(data= umapMLANAS1AllS2, aes(UMAP_1, UMAP_2, color = sampleNum), contour_var = "ndensity") +
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3DensityMemoryClusterMLANAAll.svg'), width = 6, height = 5.321)

###exampleBarcodes: toNGFR = B134; toMiddle = B72; toVCAM1 = B110; toMLANA= B81
umapMLANAS1AllS2NGFR = umapMLANAS1AllS2 %>% filter(barcodeName == "B134")
plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 1, shape = 16) + 
  geom_point(data= umapMLANAS1AllS2NGFR, aes(UMAP_1, UMAP_2, color = sampleNum), size = 2, shape = 16) + 
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3memoryS1S2MLANA_NGFR134.svg'), width = 6, height = 5.321)

umapMLANAS1AllS2Middle = umapMLANAS1AllS2 %>% filter(barcodeName == "B72")
plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 1, shape = 16) + 
  geom_point(data= umapMLANAS1AllS2Middle, aes(UMAP_1, UMAP_2, color = sampleNum), size = 2, shape = 16) + 
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3memoryS1S2MLANA_Middle72.svg'), width = 6, height = 5.321)

umapMLANAS1AllS2VCAM1 = umapMLANAS1AllS2 %>% filter(barcodeName == "B110")
plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 1, shape = 16) + 
  geom_point(data= umapMLANAS1AllS2VCAM1, aes(UMAP_1, UMAP_2, color = sampleNum), size = 2, shape = 16) + 
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3memoryS1S2MLANA_VCAM1110.svg'), width = 6, height = 5.321)

umapMLANAS1AllS2MLANA = umapMLANAS1AllS2 %>% filter(barcodeName == "B81")
plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 1, shape = 16) + 
  geom_point(data= umapMLANAS1AllS2MLANA, aes(UMAP_1, UMAP_2, color = sampleNum), size = 2, shape = 16) + 
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3memoryS1S2MLANA_MLANA81.svg'), width = 6, height = 5.321)

####
####PieChart
none = nrow(fractionColony %>% filter(sampleNum == "S1")) - length(unique(umapMLANAS1AllS2$BC50StarcodeD8))
twin = length(unique(umapMLANAS1AllS2$BC50StarcodeD8))
dataPie = data.frame(type = c("none", "twin"),
                     value = c(none/(none+twin), twin/(none+twin)))

plot = ggplot(dataPie, aes(x="", y=value, fill=type)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  theme_classic() +
  theme(axis.text.x=element_blank(), legend.position="none") +
  geom_text(aes(label = value*100), size=5)
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3pieMLANA.svg'), width = 6, height = 6)
#############################################################################
#####Where do NGFR in Trametinib go to in Veurafenib
#############################################################################
clusterMLANACellIDs = umapClusters %>% filter(scanorama_snn_res.0.5 %in% 4)
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,scanorama_snn_res.0.5)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary1 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(total = sum(total*nColony.y))
s1NGFR = colonySizeClusterMLANASummary1$total[colonySizeClusterMLANASummary1$sampleNum=="S1"]
s2NGFR = (colonySizeClusterMLANASummary1$total[colonySizeClusterMLANASummary1$sampleNum=="S3"])*(nrow(jointUMAPS1))/(nrow(jointUMAPS3))
percetEnrichment = (s2NGFR/s1NGFR)
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <6]),
                                                                                                     large = sum(total[nColony.y >5]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum)
colonySizeClusterMLANASummary2$type <- factor(colonySizeClusterMLANASummary2$type, levels=c("large", "small", "singlet"))

plot = ggplot(data=colonySizeClusterMLANASummary2, aes(x=sampleNum, y=value, fill=type)) +
  geom_bar(stat="identity") +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3NGFRSummary.svg'), width = 3, height = 5)

jointUMAPS1 = jointUMAP %>% filter(sampleNum=="S1")
jointUMAPS3 = jointUMAP %>% filter(sampleNum=="S3")

jointUMAPS1Barcodes = jointUMAP %>% filter(sampleNum=="S1") %>% select(BC50StarcodeD8) %>% unique()
jointUMAPS3Barcodes = jointUMAP %>% filter(sampleNum=="S3")  %>% select(BC50StarcodeD8) %>% unique()

jointBarcodesOnlyBoth = inner_join(jointUMAPS1Barcodes,jointUMAPS3Barcodes, by = "BC50StarcodeD8") 
jointBarcodesOnlyBoth = jointBarcodesOnlyBoth %>% mutate(barcodeName = c(1:nrow(jointBarcodesOnlyBoth)))
jointBarcodesOnlyBoth$barcodeName = sub("^", "B", jointBarcodesOnlyBoth$barcodeName)
jointBarcodesOnlyBothList = jointBarcodesOnlyBoth$BC50StarcodeD8

jointUMAPOnlyBoth = jointUMAP %>% filter(BC50StarcodeD8 %in% jointBarcodesOnlyBothList)
finalUMAPJoint = inner_join(jointUMAPOnlyBoth,jointBarcodesOnlyBoth, by = "BC50StarcodeD8") %>% select(-cellID, -nUMI)
fractionColonyMLANAS1 = fractionColony %>% filter(sampleNum == "S3")
barcodesMLANAS1 = fractionColonyMLANAS1$BC50StarcodeD8 %>% unique()

umapMLANAS1AllS3 = finalUMAPJoint %>% filter(BC50StarcodeD8 %in% barcodesMLANAS1)

plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 1, shape = 16) + 
  geom_density_2d(data= umapMLANAS1AllS3, aes(UMAP_1, UMAP_2, color = sampleNum), contour_var = "ndensity") +
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3DensityMemoryClusterNGFRAll.svg'), width = 6, height = 5.321)

plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 0.3, shape = 16) + 
  geom_point(data= umapMLANAS1AllS3, aes(UMAP_1, UMAP_2, color = sampleNum), size = 0.3, shape = 16) + 
  facet_wrap(~barcodeName, ncol = 7) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3memoryS1S2NGFRALL_individual.svg'), width = 8, height = 8)

###exampleBarcodes: toNGFR = B18; toMLANA= B145
umapMLANAS1AllS2NGFR = umapMLANAS1AllS3 %>% filter(barcodeName == "B18")
plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 1, shape = 16) + 
  geom_point(data= umapMLANAS1AllS2NGFR, aes(UMAP_1, UMAP_2, color = sampleNum), size = 2, shape = 16) + 
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3memoryS1S2NGFR_NGFR18.svg'), width = 6, height = 5.321)

umapMLANAS1AllS2MLANA = umapMLANAS1AllS3 %>% filter(barcodeName == "B145")
plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93", size = 1, shape = 16) + 
  geom_point(data= umapMLANAS1AllS2MLANA, aes(UMAP_1, UMAP_2, color = sampleNum), size = 2, shape = 16) + 
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() + theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3memoryS1S2NGFR_MLANA145.svg'), width = 6, height = 5.321)

####PieChart
none = nrow(fractionColony %>% filter(sampleNum == "S3")) - length(unique(umapMLANAS1AllS3$BC50StarcodeD8))
twin = length(unique(umapMLANAS1AllS3$BC50StarcodeD8))
dataPie = data.frame(type = c("none", "twin"),
                     value = c(none/(none+twin), twin/(none+twin)))

plot = ggplot(dataPie, aes(x="", y=value, fill=type)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  theme_classic() +
  theme(axis.text.x=element_blank(), legend.position="none") +
  geom_text(aes(label = value*100), size=5)
ggsave(plot, file = paste0(plot3Directory, 'FM02_s3s1pieNGFR.svg'), width = 6, height = 6)
#############################################################################
####Singlets and large colonies in S1 and S2 in cluster NGFR
#############################################################################

colonyS1 = fractionColony %>% filter(nColony.y >1, sampleNum == "S1") 
colonyS1UMAP = jointUMAP %>% filter(sampleNum =="S1") %>% filter(BC50StarcodeD8 %in% c(colonyS1$BC50StarcodeD8))

colonyS3 = fractionColony %>% filter(nColony.y >1, sampleNum == "S3") 
colonyS3UMAP = jointUMAP %>% filter(sampleNum =="S3") %>% filter(BC50StarcodeD8 %in% c(colonyS3$BC50StarcodeD8))

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1, shape = 16) +
  geom_point(data = colonyS1UMAP, aes(x = UMAP_1, y = UMAP_2, col=BC50StarcodeD8), size = 2, shape = 16) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3_clusterNGFRlargeColoniesS1.svg'), width = 6, height =  5.321)

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1, shape = 16) +
  geom_point(data = colonyS3UMAP, aes(x = UMAP_1, y = UMAP_2, col=BC50StarcodeD8), size = 2, shape = 16) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3_clusterNGFRlargeColoniesS3.svg'), width = 6, height =  5.321)

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93",size = 1, shape = 16) +
  geom_point(data = colonyS1UMAP, aes(x = UMAP_1, y = UMAP_2), col="hotpink3", size = 2, shape = 16) +
  geom_point(data = colonyS3UMAP, aes(x = UMAP_1, y = UMAP_2), col="turquoise3",size = 2, shape = 16) +
  create_lpr_theme()
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3_clusterNGFRlargeColoniesS1S3.svg'), width = 6, height =  5.321)

###there is no colony >1 for S2.

SingletS1 = fractionColony %>% filter(nColony.y ==1, sampleNum == "S1") 
selectS1UMAP = jointUMAP %>% filter(sampleNum =="S1") %>% filter(BC50StarcodeD8 %in% c(SingletS1$BC50StarcodeD8))

SingletS3 = fractionColony %>% filter(nColony.y ==1, sampleNum == "S3") 
selectS3UMAP = jointUMAP %>% filter(sampleNum =="S3") %>% filter(BC50StarcodeD8 %in% c(SingletS3$BC50StarcodeD8))

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1, shape = 16) +
  geom_point(data = selectS1UMAP, aes(x = UMAP_1, y = UMAP_2, col=BC50StarcodeD8), size = 2, shape = 16) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3_clusterNGFRSingletS1.svg'), width = 6, height =  5.321)

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1, shape = 16) +
  geom_point(data = selectS3UMAP, aes(x = UMAP_1, y = UMAP_2, col=BC50StarcodeD8), size = 2, shape = 16) +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3_clusterNGFRSingletS3.svg'), width = 6, height =  5.321)

plot = ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93",size = 1, shape = 16) +
  geom_point(data = selectS1UMAP, aes(x = UMAP_1, y = UMAP_2), col="hotpink3", size = 2, shape = 16) +
  geom_point(data = selectS3UMAP, aes(x = UMAP_1, y = UMAP_2), col="turquoise3",size = 2, shape = 16) +
  create_lpr_theme()
ggsave(plot, file = paste0(plot3Directory, 'FM02_s1s3_clusterNGFRSingletS1S3.svg'), width = 6, height =  5.321)