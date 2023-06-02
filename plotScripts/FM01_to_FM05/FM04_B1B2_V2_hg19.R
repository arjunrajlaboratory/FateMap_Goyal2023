#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2021_FM04/jointAnalysis/seurat/B1B2_V1_hg19.R
##############################################

library(Seurat)
library(ggplot2)
library(tidyverse)
library(fuzzyjoin)
options(future.globals.maxSize = 4000 * 1024^2)

#########
# workflow adapted from https://satijalab.org/seurat/v3.0/integration.html
#########
home1Directory <-"/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM04/"
home2Directory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM04/'
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
s1s2ScTransform <- readRDS(paste0(home1Directory,"s1s2ScTransform_50pcs_filter.rds"))
s1s2ScTransform <- FindNeighbors(object=s1s2ScTransform, dims=1:50, verbose = FALSE)
s1s2ScTransform <- FindClusters(object=s1s2ScTransform, resolution = 0.5, verbose = FALSE)

clusterUMAP = DimPlot(s1s2ScTransform, label = TRUE, label.size = 8, group.by = "seurat_clusters") + NoLegend()
sampleUMAP = DimPlot(s1s2ScTransform, label = FALSE, group.by = "orig.ident") + NoLegend()
ggsave(sampleUMAP, file = paste0(plotDirectory, 'FM04_s1s2_sampleUMAP_hg19.svg'), width = 6, height = 5.321)
ggsave(clusterUMAP, file = paste0(plotDirectory, 'FM04_s1s2_clusterUMAP_hg19.svg'), width = 6, height = 5.321)

logNormalizedCounts = as_tibble(read.table(file = paste0(home1Directory, "logNormalizedSCTCountsS1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinates = as_tibble(read.table(file = paste0(home1Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))


################3Figure 3E,F################################################
#####Just to plot the three with the right dimensions


plot <- ggplot(umapCoordinates, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes_string(color = logNormalizedCounts$NGFR), size = 1, shape = 16) +
  scale_color_gradient(low = "gray93", high = "darkblue") +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Scaled\n Log UMI counts")
ggsave(plot, file = paste0(plot3Directory, 'FM04_MLANA.svg'), width = 6, height = 6.429)

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

umiCut = 6 # Minimum UMI cutoff for reliable analysis. For the past runs I used between 3-5 for ~500k reads/sample. Here we have ~1.5mil/sample
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

length(unique(jointUMAP$BC50StarcodeD8))

#############################################################################
##############Analysis for large colonies ########################
#############################################################################

presentInBoth = inner_join(as_tibble(unique(jointUMAPS1$BC50StarcodeD8)), as_tibble(unique(jointUMAPS2$BC50StarcodeD8)))

colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizelarge = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >20)
colonySizemedium = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >9) %>% filter(nColony <21)
colonySizelargest = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >850)

colonySizelargestUMAP = inner_join(jointUMAP, colonySizelargest, by = c("BC50StarcodeD8", "sampleNum"))
colonySizelargeUMAP = inner_join(jointUMAP, colonySizelarge, by = c("BC50StarcodeD8", "sampleNum"))
colonySizemediumUMAP = inner_join(jointUMAP, colonySizemedium, by = c("BC50StarcodeD8", "sampleNum"))

colonyLarge =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonySizelargeUMAP, aes(UMAP_1, UMAP_2, color = sampleNum))+ 
  facet_wrap(~BC50StarcodeD8, ncol = 4) +
  theme_classic() + theme(legend.position = "none") 

colonyMedium =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonySizemediumUMAP, aes(UMAP_1, UMAP_2, color = sampleNum))+ 
  facet_wrap(~BC50StarcodeD8, ncol = 4) +
  theme_classic() + theme(legend.position = "none") 

######Each Sample
####Sample 1
colonySizeAllSample1 = colonySizeAll %>% filter(sampleNum == "S1")

colonySizelargeSample1 = colonySizeAllSample1 %>% filter(nColony >5)

barcodesNamed = as_tibble(colonySizelargeSample1$BC50StarcodeD8) %>% mutate(barcodeName = c(1:nrow(colonySizelargeSample1))) %>% rename(BC50StarcodeD8 = value)
barcodesNamed$barcodeName = sub("^", "B", barcodesNamed$barcodeName)

colonySizelargeSample1 = inner_join(colonySizelargeSample1, barcodesNamed, by= "BC50StarcodeD8")

colonysample1_list1 = colonySizelargeSample1 %>% filter(nColony <12)
colonysample1_list1 = inner_join(colonysample1_list1, jointUMAPS1, by = c("BC50StarcodeD8","sampleNum"))
colonysample1_list2 = colonySizelargeSample1 %>% filter(nColony >11)
colonysample1_list2 = inner_join(colonysample1_list2, jointUMAPS1, by = c("BC50StarcodeD8","sampleNum"))

colonyList1 =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonysample1_list1, aes(UMAP_1, UMAP_2))+ 
  facet_wrap(~barcodeName, ncol = 5) +
  theme_classic() + theme(legend.position = "none") 

####Ones that cover options (list 1): B1**, B26, B29**, B43, B52**, 


colonyList2 =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonysample1_list2, aes(UMAP_1, UMAP_2))+ 
  facet_wrap(~barcodeName, ncol = 5) +
  theme_classic() + theme(legend.position = "none") 
####Ones that cover options (list 2): B21, B7**, B40**, B44**

#################FOR PAPER SUPP FIGURE 1
colonySample1Paper = colonySizelargeSample1 %>% filter(barcodeName %in% c("B1", "B52", "B7", "B40", "B44"))
colonySample1Paper_plot = inner_join(colonySample1Paper, jointUMAPS1, by = c("BC50StarcodeD8","sampleNum"))

colonyPaper_Sample1 =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonySample1Paper_plot, aes(UMAP_1, UMAP_2))+ 
  facet_wrap(~barcodeName, ncol = 3) +
  theme_classic() + theme(legend.position = "none") 
ggsave(colonyPaper_Sample1, file = paste0(plotDirectory, 'FM04_colonyPaper_Sample1.svg'), width = 4, height = 3)



####Sample 2
colonySizeAllSample2 = colonySizeAll %>% filter(sampleNum == "S2")

colonySizelargeSample2 = colonySizeAllSample2 %>% filter(nColony >5)

barcodesNamed = as_tibble(colonySizelargeSample2$BC50StarcodeD8) %>% mutate(barcodeName = c(1:nrow(colonySizelargeSample2))) %>% rename(BC50StarcodeD8 = value)
barcodesNamed$barcodeName = sub("^", "B", barcodesNamed$barcodeName)

colonySizelargeSample2 = inner_join(colonySizelargeSample2, barcodesNamed, by= "BC50StarcodeD8")

colonysample2_list1 = colonySizelargeSample2 %>% filter(nColony <12)
colonysample2_list1 = inner_join(colonysample2_list1, jointUMAPS2, by = c("BC50StarcodeD8","sampleNum"))
colonysample2_list2 = colonySizelargeSample2 %>% filter(nColony >11)
colonysample2_list2 = inner_join(colonysample2_list2, jointUMAPS2, by = c("BC50StarcodeD8","sampleNum"))

colonySample2List1 =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonysample2_list1, aes(UMAP_1, UMAP_2))+ 
  facet_wrap(~barcodeName, ncol = 5) +
  theme_classic() + theme(legend.position = "none") 

####Ones that cover options (list 1): B6***,B5**, B8, B17  


colonySample2List2 =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonysample2_list2, aes(UMAP_1, UMAP_2))+ 
  facet_wrap(~barcodeName, ncol = 5) +
  theme_classic() + theme(legend.position = "none") 
####Ones that cover options (list 2): B36, B27, B30

#################FOR PAPER SUPP FIGURE 2
colonySample2Paper = colonySizelargeSample2 %>% filter(barcodeName %in% c("B6", "B5", "B36", "B27"))
colonySample2Paper_plot = inner_join(colonySample2Paper, jointUMAPS2, by = c("BC50StarcodeD8","sampleNum"))

colonyPaper_Sample2 =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonySample2Paper_plot, aes(UMAP_1, UMAP_2))+ 
  facet_wrap(~barcodeName, ncol = 3) +
  theme_classic() + theme(legend.position = "none") 
ggsave(colonyPaper_Sample2, file = paste0(plotDirectory, 'FM04_colonyPaper_Sample2.svg'), width = 4, height = 3)

#################FOR PAPER SUPP FIGURE 3
#####TWIN ONLY COLONIES
colonySizeS1 = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(sampleNum == "S1", nColony>2) 
colonySizeS2 = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(sampleNum == "S2", nColony>2) 

colonyBoth = inner_join(as_tibble(unique(colonySizeS1$BC50StarcodeD8)), as_tibble(unique(colonySizeS2$BC50StarcodeD8)))

colonyUMAPBoth = jointUMAP %>% filter(BC50StarcodeD8 %in% colonyBoth$value)  

colonyPaper_Both =ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= colonyUMAPBoth, aes(UMAP_1, UMAP_2, color = sampleNum))+ 
  facet_wrap(~BC50StarcodeD8, ncol = 3) +
  theme_classic() + theme(legend.position = "none") 
ggsave(colonyPaper_Both, file = paste0(plotDirectory, 'FM04_ccolonyPaper_Both.svg'), width = 4, height = 3)


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
coloniesToAnalyze2 = jointUMAPS2 %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >4)
colonyUMAP = inner_join(coloniesToAnalyze2,jointUMAPS2, by = c("BC50StarcodeD8")) %>% select(BC50StarcodeD8,cellID,nColony)
colonyCluster = inner_join(colonyUMAP,umapClustersS1, by = c("cellID")) %>% select(-cellID,-sampleNum,nColony) %>% ungroup()

fractionColonies = colonyCluster %>% group_by(BC50StarcodeD8,seurat_clusters) %>% summarise(nCount = length(seurat_clusters))
fractionColonies1 <- fractionColonies %>% group_by(BC50StarcodeD8) %>% filter(nCount == max(nCount))%>% slice(1)
finalFractionColonies = inner_join(fractionColonies1,coloniesToAnalyze2, by = "BC50StarcodeD8") %>% mutate(nfractionFinal = nCount/nColony) %>% select(nfractionFinal) %>% mutate(type = "colonies")

fractionColoniesSeuratClusters = fractionColonies1 %>% select(BC50StarcodeD8,seurat_clusters)
fractionColoniesSeuratClustersColony = inner_join(fractionColoniesSeuratClusters,coloniesToAnalyze2, by = "BC50StarcodeD8")
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
write.table(fractionFinalClusters, file=paste0(home1Directory,'FM04_fractionFinalClusters_S1_snn05.tsv'), col.names = TRUE, sep='\t')

summaryfractionFinalClusters = fractionFinalClusters %>% group_by(type) %>% summarise(mean = mean(nfractionFinal))

#### for S1, random = 0.202; observed = 0.389; for S2, random = 0.238; observed = 0.413)

plot = ggplot(fractionFinalClusters, aes(x = type, y=nfractionFinal, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  geom_jitter(width = 0.1, shape = 16) +
  scale_fill_manual(values=c("grey","goldenrod3")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none") +
  ylim(0,1)
ggsave(plot, file = paste0(plotDirectory, 'FM04_ColonySpecificPlotColony_V1_withDotsFinal_S2_05.svg'), width = 4, height = 8)

##################Estimation of survival fraction in MDA-MB-231-D4

UniquelyBarcodedCells_BC18 = 86000 ###from GFP sort experiment for positive cells (starting population =200,000 and GFP positive = 43.0%)
coloniesToAnalyzeS1 = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >3) %>% filter(sampleNum=="S1")
coloniesToAnalyzeS2 = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >3) %>% filter(sampleNum=="S2")

survivingClonesSample1Fraction = 1/(nrow(coloniesToAnalyzeS1)/UniquelyBarcodedCells_BC18)
survivingClonesSample2Fraction = 1/(nrow(coloniesToAnalyzeS2)/UniquelyBarcodedCells_BC18)

