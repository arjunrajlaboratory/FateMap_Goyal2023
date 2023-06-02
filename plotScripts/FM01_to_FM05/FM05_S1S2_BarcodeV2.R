#####################NOTE#####################
#Cleaned up version of V2: /Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2021_FM05_06/jointAnalysis/seurat/FM05_S1S2_Barcode.R
##############################################

library(Seurat)
library(ggplot2)
library(tidyverse)
library(fuzzyjoin)
options(future.globals.maxSize = 4000 * 1024^2)

#########
# workflow adapted from https://satijalab.org/seurat/v3.0/integration.html
#########
home1Directory <-"/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM05/"
home2Directory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM05/'
plotDirectory <- "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/plotData/analysisPlots/"

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
#############################################################################
#DON'T NEED TO RUN BLOCK THIS AGAIN
#############################################################################

s1s2ScTransform <- readRDS(paste0(home1Directory,"s1s2250nMWM983B_ScTransform_50pcs_filter.rds"))

clusterUMAP = DimPlot(s1s2ScTransform, label = TRUE, label.size = 8, group.by = "seurat_clusters") + NoLegend() 
sampleUMAP = DimPlot(s1s2ScTransform, label = FALSE, group.by = "orig.ident") + NoLegend()
ggsave(sampleUMAP, file = paste0(plotDirectory, 'FM05_s1s2_sampleUMAP_hg19.svg'), width = 6, height = 5.321)
ggsave(clusterUMAP, file = paste0(plotDirectory, 'FM05_s1s2_clusterUMAP_hg19.svg'), width = 6, height = 5.321)

#############################################################################
#DON'T NEED TO RUN BLOCK THIS AGAIN
#############################################################################

logNormalizedCounts = as_tibble(read.table(file = paste0(home1Directory, "logNormalizedSCTCountsS1S2_FM05_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinates = as_tibble(read.table(file = paste0(home1Directory, "umapCoordinatesSCT_S1S2_FM05_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))


################3Figure 3E,F################################################
#####Just to plot the three with the right dimensions


plot <- ggplot(umapCoordinates, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes_string(color = logNormalizedCounts$WNT5A), size = 1, shape = 16) +
  scale_color_gradient(low = "gray93", high = "darkblue") +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Scaled\n Log UMI counts")
ggsave(plot, file = paste0(plotDirectory, 'FM05_WNT5A.svg'), width = 6, height = 6.429)

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

umiCut = 3 # Minimum UMI cutoff for reliable analysis. For the past runs I used between 3-5 for ~500k reads/sample. Here we have ~1.5mil/sample
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

umapCoordinates_S1 = umapCoordinates %>% filter(sampleNum == "S1")
umapCoordinates_S2 = umapCoordinates %>% filter(sampleNum == "S2")

jointUMAP = inner_join(linCountTooverlaps, umapCoordinates, by = c("cellID", "sampleNum")) %>% select(-nLineages)


jointUMAPS1 = jointUMAP %>% filter(sampleNum == "S1")
jointUMAPS2 = jointUMAP %>% filter(sampleNum == "S2")

length(unique(jointUMAP$BC50StarcodeD8))
length(unique(jointUMAPS1$BC50StarcodeD8))
length(unique(jointUMAPS2$BC50StarcodeD8))

#############################################################################
##############Analysis for large colonies ########################
#############################################################################

#################FOR PAPER SUPP FIGURE 3
#####TWIN ONLY COLONIES
colonySizeS1 = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(sampleNum == "S1", nColony>1) 
colonySizeS2 = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(sampleNum == "S2", nColony>1) 
colonyBoth = inner_join(as_tibble(unique(colonySizeS1$BC50StarcodeD8)), as_tibble(unique(colonySizeS2$BC50StarcodeD8)))
barcodesNamed = as_tibble(colonyBoth$value) %>% mutate(barcodeName = c(1:nrow(colonyBoth))) %>% rename(BC50StarcodeD8 = value)
barcodesNamed$barcodeName = sub("^", "B", barcodesNamed$barcodeName)

colonyUMAPBoth =inner_join(jointUMAP, barcodesNamed, by = "BC50StarcodeD8")

colonyPaper_Both =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonyUMAPBoth, aes(UMAP_1, UMAP_2, color = sampleNum))+ 
  facet_wrap(~barcodeName, ncol = 3) +
  theme_classic() + theme(legend.position = "none") 
ggsave(colonyPaper_Both, file = paste0(plotDirectory, 'FM05_ccolonyPaperWM983B_Both.svg'), width = 4, height = 9)

###selected representation FOR PAPER: (COLONY) B7, B12, B6, B16, B1
barcodesForPaper = c("B7", "B12", "B6", "B16", "B1")
colonyUMAPBothPaper = colonyUMAPBoth %>% filter(barcodeName %in% barcodesForPaper)
colonyPaper_BothPaper =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonyUMAPBothPaper, aes(UMAP_1, UMAP_2, color = sampleNum))+ 
  facet_wrap(~barcodeName, ncol = 3) +
  theme_classic() + theme(legend.position = "none") 
ggsave(colonyPaper_BothPaper, file = paste0(plotDirectory, 'FM05_ccolonyPaperWM983B_PAPER.svg'), width = 4, height = 3)

###############
singltetSizeS1 = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(sampleNum == "S1", nColony==1) 
singltetSizeS2 = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(sampleNum == "S2", nColony==1) 
singletBoth = inner_join(as_tibble(unique(singltetSizeS1$BC50StarcodeD8)), as_tibble(unique(singltetSizeS2$BC50StarcodeD8)))
barcodesNamed = as_tibble(singletBoth$value) %>% mutate(barcodeName = c(1:nrow(singletBoth))) %>% rename(BC50StarcodeD8 = value)
barcodesNamed$barcodeName = sub("^", "B", barcodesNamed$barcodeName)

B = as_tibble(21:40)
B$value = sub("^", "B", B$value)
  
singletUMAPBoth =inner_join(jointUMAP, barcodesNamed, by = "BC50StarcodeD8")
singletUMAPBoth_Select = singletUMAPBoth %>% filter(barcodeName %in% B$value)

singletUMAPBoth_plot = ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= singletUMAPBoth_Select, aes(UMAP_1, UMAP_2, color = sampleNum))+ 
  facet_wrap(~barcodeName, ncol = 3) +
  theme_classic() + theme(legend.position = "none") 


###selected representation: (SINGLET): B39, B33, B34, B5
barcodesForPapersinglet = c("B5", "B39", "B33", "B34")
singletUMAPBothPaper = singletUMAPBoth %>% filter(barcodeName %in% barcodesForPapersinglet)
singletUMAPBothPaper_Paper =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= singletUMAPBothPaper, aes(UMAP_1, UMAP_2, color = sampleNum))+ 
  facet_wrap(~barcodeName, ncol = 3) +
  theme_classic() + theme(legend.position = "none") 
ggsave(singletUMAPBothPaper_Paper, file = paste0(plotDirectory, 'FM05_singletPaperWM983B_PAPER.svg'), width = 4, height = 3)

###########################################################################
##############################INDIVIDUAL COLONIES ##############################
##########################################################################################

###########SAMPLE 1 ########################
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAllSample1 = colonySizeAll %>% filter(sampleNum == "S1")
barcodesNamed = as_tibble(colonySizeAllSample1$BC50StarcodeD8) %>% mutate(barcodeName = c(1:nrow(colonySizeAllSample1))) %>% rename(BC50StarcodeD8 = value)
barcodesNamed$barcodeName = sub("^", "B", barcodesNamed$barcodeName)
colonySizeAllRenamedSample1 = inner_join(colonySizeAllSample1, barcodesNamed, by= "BC50StarcodeD8")

colonySizelargeSample1 = colonySizeAllRenamedSample1 %>% filter(nColony >10)
colonySizelargeSample1 = colonySizeAllRenamedSample1 %>% filter(nColony >2)
colonySizeMediumSample1 = colonySizeAllRenamedSample1 %>% filter(nColony  %in% 4:10)
colonySizeSmallSample1 = colonySizeAllRenamedSample1 %>% filter(nColony  %in% 3)
colonySizeSingletSample1 = colonySizeAllRenamedSample1 %>% filter(nColony == 1)

colonySizelargeSample1UMAP = inner_join(jointUMAP, colonySizelargeSample1, by = c("BC50StarcodeD8", "sampleNum"))
colonySizemediumSample1UMAP = inner_join(jointUMAP, colonySizeMediumSample1, by = c("BC50StarcodeD8", "sampleNum"))
colonySizemsmallSample1UMAP = inner_join(jointUMAP, colonySizeSmallSample1, by = c("BC50StarcodeD8", "sampleNum"))

colonySizesingletSample1UMAP = inner_join(jointUMAP, colonySizeSingletSample1, by = c("BC50StarcodeD8", "sampleNum"))

ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonySizesingletSample1UMAP, aes(UMAP_1, UMAP_2, color = sampleNum))+ 
  facet_wrap(~barcodeName, ncol = 4) +
  theme_classic() + theme(legend.position = "none") 

###select: selectedPaper = c("B1740", "B2382", "B4318", "B1344", "B2953", "B185", "B2267", "B2613", "B4165", "B431")
selectedPaperSample1 = c("B1740", "B2382", "B4318", "B2267", "B2613", "B431")
colonySample1Paper = colonySizeAllRenamedSample1 %>% filter(barcodeName %in% selectedPaperSample1)
colonyUMAPSample1Paper = inner_join(jointUMAP, colonySample1Paper, by = c("BC50StarcodeD8", "sampleNum"))

colonyPaper_sample1Paper =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonyUMAPSample1Paper, aes(UMAP_1, UMAP_2, color = sampleNum))+ 
  facet_wrap(~barcodeName, ncol = 3) +
  theme_classic() + theme(legend.position = "none") 
ggsave(colonyPaper_sample1Paper, file = paste0(plotDirectory, 'FM05_ccolonyPaperSAMPLE1WM983B_PAPER.svg'), width = 4, height = 3)


####All Singlets
singletPaperAll =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonySizesingletSample1UMAP, aes(UMAP_1, UMAP_2, color = sampleNum))+ 
  theme_classic() + theme(legend.position = "none") 
ggsave(singletPaperAll, file = paste0(plotDirectory, 'FM05_SingletAllSample1WM983B_PAPER.svg'), width = 6, height = 5.321)

largePaperAll =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonySizelargeSample1UMAP, aes(UMAP_1, UMAP_2, color = sampleNum))+ 
  theme_classic() + theme(legend.position = "none") 
ggsave(largePaperAll, file = paste0(plotDirectory, 'FM05_ColonyAllSample1WM983B_PAPER.svg'), width = 6, height = 5.321)

#######################Lineage Constrain WM983B
####Dmoninant Cluster analysis
s1s2ScTransform <- FindClusters(object=s1s2ScTransform, resolution = 0.5, verbose = FALSE)
jointUMAPS1 = jointUMAP %>% filter(sampleNum == "S1")
jointUMAPS2 = jointUMAP %>% filter(sampleNum == "S2")

umapClusters = (s1s2ScTransform[['seurat_clusters']])
cells_Clusters = rownames(umapClusters) #CellIds with Sample number as prefix
cells_Clusters_cellID = sub("S\\d_", "", cells_Clusters)
cells_Clusters_sampleNum = gsub("[^S12]", "", cells_Clusters)
umapClusters = as_tibble(umapClusters)
umapClusters = umapClusters %>% mutate(cellID = cells_Clusters_cellID,
                                       sampleNum = cells_Clusters_sampleNum)


umapClustersS1 = umapClusters %>% filter(sampleNum == "S2")
coloniesToAnalyze = jointUMAPS2 %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >4)
colonyUMAP = inner_join(coloniesToAnalyze,jointUMAPS2, by = c("BC50StarcodeD8")) %>% select(BC50StarcodeD8,cellID,nColony)
colonyCluster = inner_join(colonyUMAP,umapClustersS1, by = c("cellID")) %>% select(-cellID,-sampleNum,nColony)

fractionColonies = colonyCluster %>% group_by(BC50StarcodeD8,seurat_clusters) %>% summarise(nCount = length(seurat_clusters))
fractionColonies1 <- fractionColonies %>% group_by(BC50StarcodeD8) %>% filter(nCount == max(nCount))%>% slice(1)
finalFractionColonies = inner_join(fractionColonies1,coloniesToAnalyze, by = "BC50StarcodeD8") %>% mutate(nfractionFinal = nCount/nColony) %>% select(nfractionFinal) %>% mutate(type = "colonies")


fractionColoniesSeuratClusters = fractionColonies1 %>% select(BC50StarcodeD8,seurat_clusters)
fractionColoniesSeuratClustersColony = inner_join(fractionColoniesSeuratClusters,coloniesToAnalyze, by = "BC50StarcodeD8")
#fractionColoniesA$seurat_clusters <- as.numeric(fractionColoniesA$seurat_clusters) # then check that fractionColoniesA looks the same as fractionColonies
#####
maxFraction = c()
medianFraction = c()
meanFraction = c()

for (i in c(1:length(fractionColoniesSeuratClustersColony$BC50StarcodeD8))) {
  for (j in c(1:50)) {
    subSampled = sample_n(colonyCluster, fractionColoniesSeuratClustersColony$nColony[i])
    fraction = subSampled %>% group_by(seurat_clusters) %>% summarise(nFraction = length(seurat_clusters)/fractionColoniesSeuratClustersColony$nColony[i])
    if (fractionColoniesSeuratClustersColony$seurat_clusters[i] %in% fraction$seurat_clusters) {
      fractionTemp = fraction %>% filter(seurat_clusters == fractionColoniesSeuratClustersColony$seurat_clusters[i])
      maxFraction[j] = fractionTemp$nFraction
    } else {
      maxFraction[j] = 0
    }
  }
  meanFraction[i] = mean(maxFraction)
  medianFraction[i] = median(maxFraction) 
}

finalFractionRandom = tibble(nfractionFinal = numeric()) %>% add_row(nfractionFinal = meanFraction) %>% mutate(type = "random")
fractionFinalClusters = bind_rows(finalFractionColonies,finalFractionRandom)
fractionFinalClusters$type <- factor(fractionFinalClusters$type, levels=c("random", "colonies"))
write.table(fractionFinalClusters, file=paste0(home1Directory,'FM05_fractionFinalClusters_S2_snn05.tsv'), col.names = TRUE, sep='\t')

summaryfractionFinalClusters = fractionFinalClusters %>% group_by(type) %>% summarise(mean = mean(nfractionFinal))

#### for S1, random = 0.303; observed = 0.636; for S2, random = 0.372; observed = 0.578)


plot = ggplot(fractionFinalClusters, aes(x = type, y=nfractionFinal, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  geom_jitter(width = 0.1, shape = 16) +
  scale_fill_manual(values=c("grey","goldenrod3")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none") +
  ylim(0,1)
ggsave(plot, file = paste0(plotDirectory, 'FM05_ColonySpecificPlotColony_V1_withDotsFinal_S2_05.svg'), width = 4, height = 8)
