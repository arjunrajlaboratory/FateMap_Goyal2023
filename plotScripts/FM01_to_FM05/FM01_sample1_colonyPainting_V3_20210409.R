#####################NOTE#####################
#Cleaned up version of V2: /Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2019_FM01/SubmissionScripts/FINAL_COPIED/finalFiguresPaper/sample1_colonyPainting_V2_20210409.R
##############################################

library(ggplot2)
library(tidyverse)
library(reshape2)
library(svglite)
library(Seurat)

#### Memory across sampels_using data from s1s2_integration_scTransForm
home1Directory <-"/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM01/"
home2Directory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM01/'
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


#####
umapCoordinates = as_tibble(read.table(file = paste0(home1Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
logNormalizedCounts = as_tibble(read.table(file = paste0(home1Directory, "logNormalizedSCTCountsS1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
sample1_2 = readRDS(file = paste0(home1Directory,'s1s2ScTransform_50pcs_filter.rds'))  ##to read if needed.
sample1_2 <- FindNeighbors(object=sample1_2, dims=1:50, verbose = FALSE)
sample1_2 <- FindClusters(object=sample1_2, resolution = 1.2, verbose = FALSE)
saveRDS(sample1_2, file = paste0(home1Directory,'s1s2ScTransform_50pcs_filter.rds'))

sample1_2 <- FindClusters(object=sample1_2, resolution = .45, verbose = FALSE)

#######Figure 1B UMAP all ###########
clusterUMAP = DimPlot(sample1_2, reduction = "umap", label = FALSE, group.by = "integrated_snn_res.0.6") + NoLegend()
clusterUMAPWCluster = DimPlot(sample1_2, reduction = "umap", label = TRUE, group.by = "integrated_snn_res.0.6") + NoLegend()
sampleUMAP = DimPlot(sample1_2, reduction = "umap", group.by = "orig.ident")+ NoLegend()
sampleSplitUMAP = DimPlot(sample1_2, reduction = "umap",group.by = "integrated_snn_res.0.6", label = TRUE, label.size = 8, split.by = "orig.ident",) + NoLegend()

ggsave(clusterUMAP, file = paste0(plotDirectory, 'FM01_s1s2_clusterUMAP_onlyScanorma_snn06.svg'), width = 6, height = 5.25)
ggsave(clusterUMAPWCluster, file = paste0(plotDirectory, 'FM01_s1s2_clusterUMAPWNumbers_onlyScanorma_snn06.svg'), width = 6, height = 5.25)
ggsave(sampleUMAP, file = paste0(plotDirectory, 'FM01_s1s2_sampleUMAP_onlyScanorma06.svg'), width = 6, height = 5.25)
ggsave(sampleSplitUMAP, file = paste0(plotDirectory, 'FM01_s1s2_sampleSplitUMAP_onlyScanorma_snn06.svg'), width = 12, height = 5.822)

#######
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum=="S1")
logNormalizedCountsS1 = logNormalizedCounts %>% filter(sampleNum=="S1")

plot <- ggplot(umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes_string(color = logNormalizedCountsS1$CXCL3), size = 1, shape = 16) +
  scale_color_gradient(low = "gray93", high = "darkblue") +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Scaled\n Log UMI counts")
ggsave(plot, file = paste0(plotDirectory, 'FM01_UBC_s1_test.svg'), width = 6, height = 6.286)

### cacl width: (96/79.3)/(96/84.3)*6 = 6.38 ==> needs correction and settles to 6.286
####size in illustrator: 96 by 79.3

#######Barcode entries for Sample 1 and Sample 2
barcode50 = as_tibble(read.table(paste0(home2Directory, 'stepFourStarcodeShavedReads50.txt'), stringsAsFactors=F, header = T))
umiCut = 4 # Minimum UMI cutoff for reliable analysis. Can be played around with between 3 to 5.
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.  

upToLineageCounting = barcode50 %>% 
  filter(sampleNum == "1") %>%
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
#####

jointUMAP = inner_join(linCountTooverlaps, umapCoordinatesS1, by = c("cellID", "sampleNum")) %>% select(-nLineages)
singletOnlyS1 = jointUMAP %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony ==4) %>% select(-nColony)
singletOnlyS1 = jointUMAP %>% filter(BC50StarcodeD8 %in% singletOnlyS1$BC50StarcodeD8)

ggplot() + 
  geom_point(data= umapCoordinates, aes(UMAP_1, UMAP_2), color="gray93") + 
  geom_point(data= singletOnlyS1, aes(UMAP_1, UMAP_2), color="hotpink3")+ 
  facet_wrap(~BC50StarcodeD8, ncol = 4) +
  theme_classic() + theme(legend.position = "none")

####################20210318
###merging cellIDs and Barcodes
cellIDs = logNormalizedCountsS1$cellID
rowstokeep = which(cellIDs %in% unlist(linCountTooverlaps[,1])) #comparing cells present in both BCSeq and Seurat
logNormalizedCountsSubset = logNormalizedCountsS1[rowstokeep,]
logNormalizedCountsSubsetWBarcodes = inner_join(logNormalizedCountsSubset,linCountTooverlaps, by = c("cellID", "sampleNum"))


cutoff_ratio = 0.6
upperCutoff = 200
lowerCutoff = 50

ggplot(logNormalizedCountsSubsetWBarcodes, aes(x=SOX10)) +
  geom_density(alpha=0.5)

###For Figure 1, paper
#geneSet = c("ACTA2", "ACTG2", "MYOCD","TAGLN")
geneSet = c("ACTA2", "ACTG2", "MYOCD") 
geneSet = c("ACTA2", "MYOCD", "NGFR", "S100B", "MLANA", "SOX10") 
ColoniesSample1 = logNormalizedCountsSubsetWBarcodes %>% filter (sampleNum == "S1") %>% group_by(BC50StarcodeD8) %>% select(BC50StarcodeD8, ACTA2) %>% summarise(colonySize = length(BC50StarcodeD8), sum = sum(ACTA2 > 0.5)) %>% mutate(ratio = sum / colonySize) %>% filter(ratio > cutoff_ratio, colonySize < upperCutoff & colonySize > lowerCutoff) %>% mutate(gene = 'ACTA2')
BCsIndicesToKeep = which(logNormalizedCountsSubsetWBarcodes$BC50StarcodeD8 %in% unlist(ColoniesSample1[,1]))
BCsToKeep = logNormalizedCountsSubsetWBarcodes[BCsIndicesToKeep,]  %>% select(geneSet, cellID)
colonyS1Coordinates = inner_join(BCsToKeep, umapCoordinates, by = c("cellID"))

allCells = logNormalizedCountsSubsetWBarcodes %>% select(geneSet) %>% mutate(type = "all")
colony = colonyS1Coordinates %>% select(geneSet) %>% mutate(type = "l1")
densityPlot = bind_rows(allCells,colony)

meltData = melt(data = densityPlot, id.vars = c("type"), measure.vars = geneSet)
meltDataL1 = meltData

######UMAP PLOT
plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = colonyS1Coordinates, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_ACTA2Colony.svg'), width = 6, height = 5.325)


#####VCAM1 colony
cellIDs = logNormalizedCountsS1$cellID
rowstokeep = which(cellIDs %in% unlist(linCountTooverlaps[,1])) #comparing cells present in both BCSeq and Seurat
logNormalizedCountsSubset = logNormalizedCountsS1[rowstokeep,]
logNormalizedCountsSubsetWBarcodes = inner_join(logNormalizedCountsSubset,linCountTooverlaps, by = c("cellID", "sampleNum"))

ggplot(logNormalizedCountsSubsetWBarcodes, aes(x=VCAM1)) +
  geom_density(alpha=0.5)
cutoff_ratio = 0.35
upperCutoff = 200
lowerCutoff = 70

###colony 1: 
geneSet = c("VCAM1", "APOE", "FOXF2") 

ColoniesSample1 = logNormalizedCountsSubsetWBarcodes %>% filter (sampleNum == "S1") %>% group_by(BC50StarcodeD8) %>% select(BC50StarcodeD8, VCAM1) %>% summarise(colonySize = length(BC50StarcodeD8), sum = sum(VCAM1 > 0.5)) %>% mutate(ratio = sum / colonySize) %>% filter(ratio > cutoff_ratio, colonySize < upperCutoff & colonySize > lowerCutoff) %>% mutate(gene = 'VCAM1')
####colonies: 1: "ATTCGTGGACGACTTCAACTTCAACGTGATGGACGAGCAGTTGCTGATGT" ; 2: "ATTCGAGTTCCTGTTGTAGGTCTTGCAGTAGTTCCAGGACTAGCACAAGG"; 3: "ATTGGACCAGTTGATCAAGTACGTCATGTTCTTGGACGACGTCGTCCTCA"

ColoniesSample1 = ColoniesSample1 %>% filter(BC50StarcodeD8 == "ATTCGTGGACGACTTCAACTTCAACGTGATGGACGAGCAGTTGCTGATGT")
ColoniesSample1 = ColoniesSample1 %>% filter(BC50StarcodeD8 == "ATTCGAGTTCCTGTTGTAGGTCTTGCAGTAGTTCCAGGACTAGCACAAGG") 
ColoniesSample1 = ColoniesSample1 %>% filter(BC50StarcodeD8 == "ATTGGACCAGTTGATCAAGTACGTCATGTTCTTGGACGACGTCGTCCTCA") 

BCsIndicesToKeep = which(logNormalizedCountsSubsetWBarcodes$BC50StarcodeD8 %in% unlist(ColoniesSample1[,1]))
BCsToKeep = logNormalizedCountsSubsetWBarcodes[BCsIndicesToKeep,]  %>% select(geneSet, cellID)
colonyS1Coordinates = inner_join(BCsToKeep, umapCoordinates, by = c("cellID"))

######UMAP PLOT
plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = colonyS1Coordinates, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_VCAM1Colony3.svg'), width = 6, height = 5.325)
##################

######NGFR
cutoff_ratio = 0.6
upperCutoff = 15
lowerCutoff = 9

cutoff_ratio = 0.6
upperCutoff = 200
lowerCutoff = 50

ggplot(logNormalizedCountsSubsetWBarcodes, aes(x=AXL)) +
  geom_density(alpha=0.5)

geneSet = c("NGFR", "S100B", "ITGA6") 
geneSet = c("ACTA2", "MYOCD", "NGFR", "S100B", "MLANA", "SOX10") 
ColoniesSample1 = logNormalizedCountsSubsetWBarcodes %>% filter (sampleNum == "S1") %>% group_by(BC50StarcodeD8) %>% select(BC50StarcodeD8, NGFR) %>% summarise(colonySize = length(BC50StarcodeD8), sum = sum(NGFR > 0.5)) %>% mutate(ratio = sum / colonySize) %>% filter(ratio > cutoff_ratio, colonySize < upperCutoff & colonySize > lowerCutoff) %>% mutate(gene = 'NGFR')
BCsIndicesToKeep = which(logNormalizedCountsSubsetWBarcodes$BC50StarcodeD8 %in% unlist(ColoniesSample1[,1]))
BCsToKeep = logNormalizedCountsSubsetWBarcodes[BCsIndicesToKeep,]  %>% select(geneSet, cellID)
colonyS1Coordinates = inner_join(BCsToKeep, umapCoordinates, by = c("cellID"))

allCells = logNormalizedCountsSubsetWBarcodes %>% select(geneSet) %>% mutate(type = "all")
colony = colonyS1Coordinates %>% select(geneSet) %>% mutate(type = "l2")
densityPlot = bind_rows(allCells,colony)

meltData = melt(data = densityPlot, id.vars = c("type"), measure.vars = geneSet)
meltDataL2 = meltData


######UMAP PLOT
plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = colonyS1Coordinates, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_NGFRColony.svg'), width = 6, height = 5.325)

######AXL
cutoff_ratio = 0.5
upperCutoff = 60
lowerCutoff = 10

ggplot(logNormalizedCountsSubsetWBarcodes, aes(x=AXL)) +
  geom_density(alpha=0.5)

geneSet = c("AXL", "SERPINE1") 
ColoniesSample1 = logNormalizedCountsSubsetWBarcodes %>% filter (sampleNum == "S1") %>% group_by(BC50StarcodeD8) %>% select(BC50StarcodeD8, AXL) %>% summarise(colonySize = length(BC50StarcodeD8), sum = sum(AXL > 0.5)) %>% mutate(ratio = sum / colonySize) %>% filter(ratio > cutoff_ratio, colonySize < upperCutoff & colonySize > lowerCutoff) %>% mutate(gene = 'AXL')
barcode1 = "ATACGTCGTGAACTTCATCTTGGTGTTCTAGTTCCTCTACTACGACATGG"
barcode2 = "ATAGTTCCAGCACCTGGAGTAGGAGAACCTGTTCATCGTGCTGGAGTTGA"
barcode3 = "ATTCCTGGAGCTCCTGTTGTAGAAGCTGTACCTGTAGGTGTTGCTGTTGG"
jointUMAPBarcode1 = jointUMAP %>% filter(BC50StarcodeD8 %in% barcode1)
jointUMAPBarcode2 = jointUMAP %>% filter(BC50StarcodeD8 %in% barcode2)
jointUMAPBarcode3 = jointUMAP %>% filter(BC50StarcodeD8 %in% barcode3)




######UMAP PLOT
plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = jointUMAPBarcode1, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_AXLColonyBarcode1.svg'), width = 6, height = 5.325)

######SOX10

cutoff_ratio = 0.6
upperCutoff = 2
lowerCutoff = 0
geneSet = c("SOX10", "MLANA", "MITF") 
geneSet = c("ACTA2", "MYOCD", "NGFR", "S100B", "MLANA", "SOX10") 

ColoniesSample1 = logNormalizedCountsSubsetWBarcodes %>% filter (sampleNum == "S1") %>% group_by(BC50StarcodeD8) %>% select(BC50StarcodeD8, SOX10) %>% summarise(colonySize = length(BC50StarcodeD8), sum = sum(SOX10 > 0.5)) %>% mutate(ratio = sum / colonySize) %>% filter(ratio > cutoff_ratio, colonySize < upperCutoff & colonySize > lowerCutoff) %>% mutate(gene = 'SOX10')
BCsIndicesToKeep = which(logNormalizedCountsSubsetWBarcodes$BC50StarcodeD8 %in% unlist(ColoniesSample1[,1]))
BCsToKeep = logNormalizedCountsSubsetWBarcodes[BCsIndicesToKeep,]  %>% select(geneSet, cellID)
colonyS1Coordinates = inner_join(BCsToKeep, umapCoordinates, by = c("cellID"))
allCells = logNormalizedCountsSubsetWBarcodes %>% select(geneSet) %>% mutate(type = "all")
colony = colonyS1Coordinates %>% select(geneSet) %>% mutate(type = "l3")
densityPlot = bind_rows(allCells,colony)
meltData = melt(data = densityPlot, id.vars = c("type"), measure.vars = geneSet)
meltDataL3 = meltData

#selectedSingletLineages = sample_n(colonyS1Coordinates, 3)

selectedSingletLineages = colonyS1Coordinates %>% filter(cellID %in% c("CGGAATTAGAGGCGTT", "TGTTGGATCTCCGCAT","TCATGCCAGGGCAGAG"))

allCells = logNormalizedCountsSubsetWBarcodes %>% select(geneSet) %>% mutate(type = "all")
colony = selectedSingletLineages %>% select(geneSet) %>% mutate(type = "selected")
densityPlot = bind_rows(allCells,colony)  
meltData = melt(data = densityPlot, id.vars = c("type"), measure.vars = geneSet)
meltDataL4 = meltData

plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = colonyS1Coordinates, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts")
ggsave(plot, file = paste0(plotDirectory, 'FM01_SingletsAll.svg'), width = 6, height = 5.325)

plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = selectedSingletLineages, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_SingletsSample.svg'), width = 6, height = 5.325)

plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = colonyS1Coordinates, aes(x = UMAP_1, y = UMAP_2, color = cellID), size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "none",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_SingletsAllColored.svg'), width = 6, height = 5.325)




###########Dot plot Arjun Asked################
meltDataCombined = bind_rows(meltDataL1, meltDataL2, meltDataL3)
meltDataCombined = bind_rows(meltDataL1, meltDataL2, meltDataL3, meltDataL4)

set.seed(2398) # any seed; makes the jitter reproducible
plot = ggplot() +
  geom_violin(data = meltDataCombined, aes(type, value, fill = type)) +
  scale_fill_manual(values=c("gray36", "hotpink3", "hotpink3", "hotpink3")) +
  facet_wrap(~variable,nrow=6, scales = 'free_y') + 
  scale_y_continuous(breaks = round(seq(min(meltData$value), max(meltData$value), by = 2),1)) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot2Directory, 'FM01_figure1Violin.svg'), width = 8, height = 6)

plot = ggplot() +
  geom_boxplot(data = meltDataCombined, aes(type, value, fill = type)) +
  scale_fill_manual(values=c("gray27", "hotpink3", "hotpink3", "gray77","green")) +
  facet_wrap(~variable,nrow=6, scales = 'free_y') + 
  scale_y_continuous(breaks = round(seq(min(meltData$value), max(meltData$value), by = 2),1)) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plot2Directory, 'FM01_figure1Box.svg'), width = 8, height = 6)

##########
set.seed(2398) # any seed; makes the jitter reproducible ###Used for Final figure
plot = ggplot() +
  geom_jitter(data = meltDataCombined, aes(type, value, colour = type), width = 0.25, height = 0.25, shape = 16) +
  geom_boxplot(data = meltDataCombined, aes(type, value)) +
  scale_color_manual(values=c("gray27", "hotpink3", "hotpink3", "gray77","green")) +
  facet_wrap(~variable,nrow=6, scales = 'free_y') + 
  scale_y_continuous(breaks = round(seq(min(meltData$value), max(meltData$value), by = 2),1)) +
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none", strip.background = element_blank(),strip.text.x = element_blank())
ggsave(plot, file = paste0(plotDirectory, 'FM01_figure1G_jitterBoxSelected.svg'), width = 8, height = 6)

####################################################################################################################################
######################################## Global cluster Analysis for colonies #######################################################
####################################################################################################################################
sample1_2 <- FindClusters(object=sample1_2, resolution = 0.8, verbose = FALSE)

jointUMAP = inner_join(linCountTooverlaps, umapCoordinatesS1, by = c("cellID", "sampleNum")) %>% select(-nLineages)

umapClusters = (sample1_2[['seurat_clusters']])
cells_Clusters = rownames(umapClusters) #CellIds with Sample number as prefix
cells_Clusters_cellID = sub("S\\d_", "", cells_Clusters)
cells_Clusters_sampleNum = gsub("[^S12]", "", cells_Clusters)
umapClusters = as_tibble(umapClusters)
umapClusters = umapClusters %>% mutate(cellID = cells_Clusters_cellID,
                                       sampleNum = cells_Clusters_sampleNum)

write.table(umapClusters, file=paste0(home1Directory,'umapClusters_s1s2Scanorama_50pcs_filter_snn0_8.tsv'), col.names = TRUE, sep='\t')

###Colony size >=5
umapClustersS1 = umapClusters %>% filter(sampleNum == "S1")
coloniesToAnalyze = jointUMAP %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >4)
colonyUMAP = inner_join(coloniesToAnalyze,jointUMAP, by = c("BC50StarcodeD8")) %>% select(BC50StarcodeD8,cellID)
colonyCluster = inner_join(colonyUMAP,umapClustersS1, by = c("cellID")) %>% select(-cellID,-sampleNum)

samplingLengths = coloniesToAnalyze$nColony
maxFraction = c()
secondMaxFraction = c()

for (i in c(1:length(samplingLengths))) {
  subSampled = sample_n(colonyCluster, samplingLengths[i])
  fraction = subSampled %>% group_by(seurat_clusters) %>% summarise(nFraction = length(seurat_clusters)/samplingLengths[i])
  maxFraction[i] = max(fraction$nFraction)
  secondMaxFraction[i] = max(fraction$nFraction[fraction$nFraction!=max(fraction$nFraction)]) 
}

secondMaxFraction[!is.finite(secondMaxFraction)] <- NA
secondMaxFraction[which(is.na(secondMaxFraction))] = maxFraction[which(is.na(secondMaxFraction))]

secondMaxFraction<-secondMaxFraction[!is.na(secondMaxFraction)]
mean(maxFraction)
mean(secondMaxFraction)

finalFractionRandom = tibble(nfractionFinal = numeric()) %>% add_row(nfractionFinal = maxFraction) %>% mutate(type = "random")
finalFractionRandomBoth = tibble(nfractionFinal = numeric()) %>% add_row(nfractionFinal = maxFraction + secondMaxFraction) %>% mutate(type = "random")

fractionColonies = colonyCluster %>% group_by(BC50StarcodeD8,seurat_clusters) %>% summarise(nCount = length(seurat_clusters))
fractionColonies1 = fractionColonies %>% group_by(BC50StarcodeD8) %>% summarise(nfraction = max(nCount), nfractionSecond = max(nCount[nCount!=max(nCount)]))
fractionColonies1$nfractionSecond[!is.finite(fractionColonies1$nfractionSecond)] <- 0 #####correcting for second maximum which does not exist

finalFractionColonies = inner_join(fractionColonies1,coloniesToAnalyze, by = "BC50StarcodeD8") %>% mutate(nfractionFinal = nfraction/nColony) %>% select(nfractionFinal) %>% mutate(type = "colonies")
finalFractionColoniesBoth = inner_join(fractionColonies1,coloniesToAnalyze, by = "BC50StarcodeD8") %>% mutate(nfractionFinal = (nfraction+nfractionSecond)/nColony) %>% select(nfractionFinal) %>% mutate(type = "colonies")


fractionFinalClusters = bind_rows(finalFractionColonies,finalFractionRandom)
fractionFinalClustersBoth = bind_rows(finalFractionColoniesBoth,finalFractionRandomBoth)

fractionFinalClusters$type <- factor(fractionFinalClusters$type, levels=c("random", "colonies"))
fractionFinalClustersBoth$type <- factor(fractionFinalClustersBoth$type, levels=c("random", "colonies"))


plot = ggplot(fractionFinalClusters, aes(x = type, y=nfractionFinal, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  scale_fill_manual(values=c("grey","goldenrod3")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none") +
  ylim(0,1)
ggsave(plot, file = paste0(plotDirectory, 'FM01_OverallPlotColony_V1.svg'), width = 5, height = 8)
ggsave(plot, file = paste0(plotDirectory, 'FM01_OverallPlotColony_V1_snn_04.svg'), width = 5, height = 8)
ggsave(plot, file = paste0(plotDirectory, 'FM01_OverallPlotColony_V1_snn_06.svg'), width = 5, height = 8)
ggsave(plot, file = paste0(plotDirectory, 'FM01_OverallPlotColony_V1_snn_12.svg'), width = 5, height = 8)

plot = ggplot(fractionFinalClustersBoth, aes(x = type, y=nfractionFinal, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  scale_fill_manual(values=c("grey","goldenrod3")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none") +
  ylim(0,1)
ggsave(plot, file = paste0(plotDirectory, 'FM01_OverallPlotColonyTopTwoClusters_V1.svg'), width = 5, height = 8)
ggsave(plot, file = paste0(plotDirectory, 'FM01_OverallPlotColonyTopTwoClusters_V1_snn_04.svg'), width = 5, height = 8)
ggsave(plot, file = paste0(plotDirectory, 'FM01_OverallPlotColonyTopTwoClusters_V1_snn_06.svg'), width = 5, height = 8)
ggsave(plot, file = paste0(plotDirectory, 'FM01_OverallPlotColonyTopTwoClusters_V1_snn_12.svg'), width = 5, height = 8)


####################################################################################################################################
######################################## Colony Specific cluster Analysis #########################################################
####################################################################################################################################
sample1_2 <- FindClusters(object=sample1_2, resolution = 0.8, verbose = FALSE)
umapClusters = as_tibble(read.table(file = paste0(home1Directory, "umapClusters_s1s2Scanorama_50pcs_filter_snn0_8.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapClusters = (sample1_2[['seurat_clusters']])
cells_Clusters = rownames(umapClusters) #CellIds with Sample number as prefix
cells_Clusters_cellID = sub("S\\d_", "", cells_Clusters)
cells_Clusters_sampleNum = gsub("[^S12]", "", cells_Clusters)
umapClusters = as_tibble(umapClusters)
umapClusters = umapClusters %>% mutate(cellID = cells_Clusters_cellID,
                                       sampleNum = cells_Clusters_sampleNum)

umapClustersS1 = umapClusters %>% filter(sampleNum == "S1")
coloniesToAnalyze = jointUMAP %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >4)
colonyUMAP = inner_join(coloniesToAnalyze,jointUMAP, by = c("BC50StarcodeD8")) %>% select(BC50StarcodeD8,cellID,nColony)
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
write.table(fractionFinalClusters, file=paste0(home1Directory,'FM01_fractionFinalClusters_snn08.tsv'), col.names = TRUE, sep='\t')

summaryfractionFinalClusters = fractionFinalClusters %>% group_by(type) %>% summarise(mean = mean(nfractionFinal))

plot = ggplot(fractionFinalClusters, aes(x = type, y=nfractionFinal, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  geom_jitter(width = 0.1, shape = 16) +
  scale_fill_manual(values=c("grey","goldenrod3")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none") +
  ylim(0,1)
ggsave(plot, file = paste0(plotDirectory, 'FM01_ColonySpecificPlotColony_V1_withDotsFinal.svg'), width = 4, height = 8)
ggsave(plot, file = paste0(plotDirectory, 'FM01_ColonySpecificPlotColony_V1_snn_04.svg'), width = 4, height = 8)
ggsave(plot, file = paste0(plotDirectory, 'FM01_ColonySpecificPlotColony_V1_snn_06.svg'), width = 4, height = 8)
ggsave(plot, file = paste0(plotDirectory, 'FM01_ColonySpecificPlotColony_V1_snn_12.svg'), width = 4, height = 8)

######representative UMAPS##############
umapClustersS1UMAP = inner_join(umapClustersS1, umapCoordinatesS1, by = c("cellID", "sampleNum"))
finalFractionColoniesAll = inner_join(fractionColonies,coloniesToAnalyze, by = "BC50StarcodeD8") %>% mutate(nfractionFinal = nCount/nColony) %>% select(nCount,nfractionFinal,BC50StarcodeD8, seurat_clusters)

umapCoordinatesB_high = jointUMAP %>% filter(BC50StarcodeD8 == "ATAGATCGACAAGGAGATGTTCTTCGTGGACATCTAGATCAACATCAAGT")  
umapCoordinatesB_medium = jointUMAP %>% filter(BC50StarcodeD8 == "ATTCTAGTTGTAGTACTAGATGATCATCATGTTGTTCTTGGTCTTGTTCA")  
umapCoordinatesB_low = jointUMAP %>% filter(BC50StarcodeD8 == "ATTCCTGATGGTCATGGAGAACCTGCTGGTCTTCCACATGGACAAGGAGC")  

B_high =  finalFractionColoniesAll %>% filter(BC50StarcodeD8 == "ATAGATCGACAAGGAGATGTTCTTCGTGGACATCTAGATCAACATCAAGT") %>% arrange(-nCount)
B_medium = finalFractionColoniesAll %>% filter(BC50StarcodeD8 == "ATTCTAGTTGTAGTACTAGATGATCATCATGTTGTTCTTGGTCTTGTTCA") %>% arrange(-nCount)
B_low = finalFractionColoniesAll %>% filter(BC50StarcodeD8 == "ATTCCTGATGGTCATGGAGAACCTGCTGGTCTTCCACATGGACAAGGAGC") %>% arrange(-nCount)
toPlotBHigh = umapClustersS1UMAP %>% filter(seurat_clusters %in% B_high$seurat_clusters[1:2])
toPlotBMedium = umapClustersS1UMAP %>% filter(seurat_clusters %in% B_medium$seurat_clusters[1:2])
toPlotBLow = umapClustersS1UMAP %>% filter(seurat_clusters %in% B_low$seurat_clusters[1:2])

#toPlotBHigh = inner_join(B_high,jointUMAP, by = c("BC50StarcodeD8")) %>% select(-nUMI,sampleNum,BC50StarcodeD8) %>% filter(seurat_clusters %in% B_high$seurat_clusters[1:2])
#toPlotBMedium = inner_join(finalFractionColoniesAll,jointUMAP, by = c("BC50StarcodeD8")) %>% select(-nUMI,sampleNum,BC50StarcodeD8) %>% filter(seurat_clusters %in% B_medium$seurat_clusters[1:2])
#toPlotBLow = inner_join(B_low,jointUMAP, by = c("BC50StarcodeD8")) %>% select(-nUMI,sampleNum,BC50StarcodeD8) %>% filter(seurat_clusters %in% B_low$seurat_clusters[1:2])

plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = toPlotBMedium, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters), size = 1.5, shape = 16) +
  scale_color_manual(values=c("gray40", "gray60")) +
  geom_point(data = umapCoordinatesB_medium, aes(x = UMAP_1, y = UMAP_2),color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_MediumClusterFraction.svg'), width = 6, height = 6.158)

plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = toPlotBHigh, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters), size = 1.5, shape = 16) +
  scale_color_manual(values=c("gray40", "gray60")) +
  geom_point(data = umapCoordinatesB_high, aes(x = UMAP_1, y = UMAP_2),color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_HighClusterFraction.svg'), width = 6, height = 6.158)

plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = toPlotBLow, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters), size = 1.5, shape = 16) +
  scale_color_manual(values=c("gray40", "gray60")) +
  geom_point(data = umapCoordinatesB_low, aes(x = UMAP_1, y = UMAP_2),color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_LowClusterFraction.svg'), width = 6, height = 6.158)


####################################################################################################################################
######################################## Neighbor analysis for bias #########################################################
####################################################################################################################################
library(RANN)
library(tripack)
library(reshape2)

pcaCoordinates = as_tibble(read.table(file = paste0(home1Directory, "pcaCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
jointPCA = inner_join(linCountTooverlaps, pcaCoordinates, by = c("cellID", "sampleNum")) %>% select(-nLineages)
jointPCAS1 = jointPCA %>% filter(sampleNum=="S1") 
PCAtoAnalyze = jointPCAS1 %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >20) ####since for smaller colonies neighbor analysis not ideal or have to reduce the number of allowed neighbors (which is ok, maybe I can do a response curve)
jointPCAS1_onlyColonies = inner_join(jointPCAS1,PCAtoAnalyze, by = "BC50StarcodeD8")

neighborRandomFinalExtract = c();
neighborColonyFinalExtract = c();

for (i in c(1:length(PCAtoAnalyze$BC50StarcodeD8))) {
  neighborRandomFinal = c();
  neighborColonyFinal = c();
  for (j in c(1:5)) {
    subSampled1 = sample_n(jointPCAS1_onlyColonies, PCAtoAnalyze$nColony[i]) %>% mutate(name = "random1")
    subSampled2 = sample_n(jointPCAS1_onlyColonies, PCAtoAnalyze$nColony[i]) %>% mutate(name = "random2")
    random12 = bind_rows(subSampled1,subSampled2)
    random12BarcodeX = random12 %>% mutate(num = c(1:nrow(random12)))
    Random1index = random12BarcodeX %>% filter(name == "random1") %>% select(num)
    Random2index = random12BarcodeX %>% filter(name == "random2") %>% select(num)
    knnPCARandom = nn2(random12[,5:54], random12[,5:54],k = min(10, nrow(subSampled1)))
    neighborKNNRandom = as_tibble(knnPCARandom[[1]]) ####gets the neighbor
    
    barcodeColony = jointPCAS1_onlyColonies %>% filter(BC50StarcodeD8 == PCAtoAnalyze$BC50StarcodeD8[i]) %>% mutate(name = "colony")
    random1Colony = bind_rows(subSampled1,barcodeColony)
    random1ColonyBarcodeX = random1Colony %>% mutate(num = c(1:nrow(random1Colony)))
    Randomindex = random1ColonyBarcodeX %>% filter(name == "random1") %>% select(num)
    Colonyindex = random1ColonyBarcodeX %>% filter(name == "colony") %>% select(num)
    knnPCAColony = nn2(random1Colony[,5:54], random1Colony[,5:54], k = min(10, nrow(random1Colony)))
    neighborKNNColony = as_tibble(knnPCAColony[[1]]) ####gets the neighbor
    
    nRandom1 = c();
    nRandom2 = c();
    nColony1 = c();
    nColony2 = c();
    
    for (k in c(1:nrow(random12))) {
      nRandom1[k] =  sum(neighborKNNRandom[k,2:pmin(10,nrow(subSampled1))] %in% Random1index$num) 
      nRandom2[k] =  sum(neighborKNNRandom[k,2:pmin(10,nrow(subSampled1))] %in% Random2index$num)
      nColony1[k] =  sum(neighborKNNColony[k,2:pmin(10,nrow(random1Colony))] %in% Randomindex$num) 
      nColony2[k] =  sum(neighborKNNColony[k,2:pmin(10,nrow(random1Colony))] %in% Colonyindex$num) 
    }
    #neighborRandom = tibble(nRandom1,nRandom2) %>% mutate(fractionR1 = nRandom1/(nRandom1+nRandom2), fractionR2= nRandom2/(nRandom1+nRandom2)) %>% mutate(num = c(1:nrow(random12)))
    neighborRandomR1 = tibble(nRandom1,nRandom2) %>% mutate(fractionR1 = nRandom1/(nRandom1+nRandom2), fractionR2= nRandom2/(nRandom1+nRandom2)) %>% mutate(num = c(1:nrow(random12))) %>% filter(num %in% Random1index$num)
    neighborRandomR2 = tibble(nRandom1,nRandom2) %>% mutate(fractionR1 = nRandom1/(nRandom1+nRandom2), fractionR2= nRandom2/(nRandom1+nRandom2)) %>% mutate(num = c(1:nrow(random12))) %>% filter(num %in% Random2index$num)
    R1R1 = as_tibble(neighborRandomR1$fractionR1) %>% mutate(isSelf = "withSelf")
    R2R2 = as_tibble(neighborRandomR2$fractionR2) %>% mutate(isSelf = "withSelf")
    R1R2 = as_tibble(neighborRandomR1$fractionR2) %>% mutate(isSelf = "withOther")
    R2R1 = as_tibble(neighborRandomR2$fractionR1) %>% mutate(isSelf = "withOther")
    toPlotRandom = bind_rows(R1R1,R2R2,R1R2,R2R1) %>% rename(fraction = value) %>% mutate(isSelf = factor(isSelf, levels=c("withSelf","withOther")))
    withSelf = toPlotRandom %>% filter(isSelf == "withSelf") %>% select(-isSelf)
    withOther = toPlotRandom %>% filter(isSelf == "withOther") %>% select(-isSelf)
    neighborRandomFinal[j]= mean(withOther$fraction)/mean(withSelf$fraction)
    
    #neighborColony = tibble(nColony1,nColony2) %>% mutate(fractionR1 = nColony1/(nColony1+nColony2), fractionR2= nColony2/(nColony1+nColony2)) %>% mutate(num = c(1:nrow(random1Colony)))
    neighborColonyR1 = tibble(nColony1,nColony2) %>% mutate(fractionR1 = nColony1/(nColony1+nColony2), fractionR2= nColony2/(nColony1+nColony2)) %>% mutate(num = c(1:nrow(random1Colony))) %>% filter(num %in% Randomindex$num)
    neighborColonyR2 = tibble(nColony1,nColony2) %>% mutate(fractionR1 = nColony1/(nColony1+nColony2), fractionR2= nColony2/(nColony1+nColony2)) %>% mutate(num = c(1:nrow(random1Colony)))%>% filter(num %in% Colonyindex$num)
    R1R1_colony = as_tibble(neighborColonyR1$fractionR1) %>% mutate(isSelf = "withSelf")
    R2R2_colony = as_tibble(neighborColonyR2$fractionR2) %>% mutate(isSelf = "withSelf")
    R1R2_colony = as_tibble(neighborColonyR1$fractionR2) %>% mutate(isSelf = "withOther")
    R2R1_colony = as_tibble(neighborColonyR2$fractionR1) %>% mutate(isSelf = "withOther")
    toPlotColony = bind_rows(R1R1_colony,R2R2_colony,R1R2_colony,R2R1_colony) %>% rename(fraction = value) %>% mutate(isSelf = factor(isSelf, levels=c("withSelf","withOther")))
    withSelfColony = toPlotColony %>% filter(isSelf == "withSelf") %>% select(-isSelf)
    withOtherColony = toPlotColony %>% filter(isSelf == "withOther") %>% select(-isSelf)
    neighborColonyFinal[j]= mean(withOtherColony$fraction)/mean(withSelfColony$fraction)
  }
  neighborRandomFinalExtract[i] = mean(neighborRandomFinal)
  neighborColonyFinalExtract[i] = mean(neighborColonyFinal)
}

neighborRandomFinalExtract1 = tibble(mixingFinal = numeric()) %>% add_row(mixingFinal = neighborRandomFinalExtract) %>% mutate(type = "random")
neighborColonyFinalExtract1 = tibble(mixingFinal = numeric()) %>% add_row(mixingFinal =  neighborColonyFinalExtract) %>% mutate(type = "colony")

plotAll = bind_rows(neighborRandomFinalExtract1,neighborColonyFinalExtract1)
plotAll$type <- factor(plotAll$type, levels=c("random", "colony"))
plot = ggplot(plotAll, aes(x = type, y=mixingFinal, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  scale_fill_manual(values=c("grey","goldenrod3")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none") +
  ylim(0,1.2)
ggsave(plot, file = paste0(plotDirectory, 'FM01_neighborAnalysisToRandom.svg'), width = 4, height = 8)

plot = ggplot(plotAll, aes(x = type, y=mixingFinal, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  geom_jitter(width = 0.1, shape = 16) +
  scale_fill_manual(values=c("grey","goldenrod3")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none") +
  ylim(0,1.2)
ggsave(plot, file = paste0(plotDirectory, 'FM01_neighborAnalysisToRandomWDots.svg'), width = 4, height = 8)

####################################################################################################################################
######################################## ColonySizeAnalysis #########################################################
####################################################################################################################################
###MLANA: 0; NGFR = 2, 12; ACTA2: 7; VCAM1 : 5,13; 10: IFIT2 10 (snn 0.4)

################################################
############this block only once################
################################################
sample1_2 <- FindClusters(object=sample1_2, resolution = 0.6, verbose = FALSE)
jointUMAP = inner_join(linCountTooverlaps, umapCoordinatesS1, by = c("cellID", "sampleNum")) %>% select(-nLineages)
umapClusters = (sample1_2[['seurat_clusters']])
cells_Clusters = rownames(umapClusters) #CellIds with Sample number as prefix
cells_Clusters_cellID = sub("S\\d_", "", cells_Clusters)
cells_Clusters_sampleNum = gsub("[^S12]", "", cells_Clusters)
umapClusters = as_tibble(umapClusters)
umapClusters = umapClusters %>% mutate(cellID = cells_Clusters_cellID,
                                       sampleNum = cells_Clusters_sampleNum)
write.table(umapClusters, file=paste0(home1Directory,'umapClusters_s1s2Scanorama_50pcs_filter_snn0_6.tsv'), col.names = TRUE, sep='\t')
###this block ^^ only once###########

umapClusters = as_tibble(read.table(file = paste0(home1Directory, "umapClusters_s1s2Scanorama_50pcs_filter_snn0_6.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))

clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(0,3))
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <4]),
                                                                                                     large = sum(total[nColony.y >3]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum)
MLANASummary = colonySizeClusterMLANASummary2 %>% mutate(cluster = "mlana") %>% mutate(percent = 100*value/(sum(value)))

####NGFR
clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(7))
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <4]),
                                                                                                     large = sum(total[nColony.y >3]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum) %>% mutate(percent = 100*value/(sum(value)))
NGFRSummary = colonySizeClusterMLANASummary2 %>% mutate(cluster = "ngfr")

####IFIT2
clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(12))
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <4]),
                                                                                                     large = sum(total[nColony.y >3]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum) %>% mutate(percent = 100*value/(sum(value)))
IFITSummary = colonySizeClusterMLANASummary2 %>% mutate(cluster = "IFIT")

####VCAM1

clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(15,4,6))
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <4]),
                                                                                                     large = sum(total[nColony.y >3]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum) %>% mutate(percent = 100*value/(sum(value)))
VCAM1Summary = colonySizeClusterMLANASummary2 %>% mutate(cluster = "VCAM1")

####ACTA2
clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(8))
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <4]),
                                                                                                     large = sum(total[nColony.y >3]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum)
ACTA2Summary = colonySizeClusterMLANASummary2 %>% mutate(cluster = "ACTA2") %>% mutate(percent = 100*value/(sum(value)))


####AXL, SERPINE1
clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(11,1,2))
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <4]),
                                                                                                     large = sum(total[nColony.y >3]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum)
AXLSummary = colonySizeClusterMLANASummary2 %>% mutate(cluster = "AXL") %>% mutate(percent = 100*value/(sum(value)))

#########
All = bind_rows(AXLSummary, ACTA2Summary,MLANASummary, NGFRSummary,IFITSummary, VCAM1Summary)
All$type <- factor(All$type, levels=c("large", "small", "singlet"))
All$cluster <- factor(All$cluster, levels=c("mlana", "ngfr", "IFIT", "ACTA2", "AXL", "VCAM1"))

write.table(All, file=paste0(plotDirectory,'FM01_AllColoniesSummary_PercentValue.tsv'), col.names = TRUE, sep='\t')

plot = ggplot(data=All, aes(x=cluster, y=percent, fill=type)) +
  geom_bar(stat="identity") +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plotDirectory, 'FM01_AllColoniesSummary_Percent.svg'), width = 8, height = 5)

plot = ggplot(data=All, aes(x=cluster, y=value, fill=type)) +
  geom_bar(stat="identity") +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plotDirectory, 'FM01_AllColoniesSummary_Value.svg'), width = 8, height = 5)