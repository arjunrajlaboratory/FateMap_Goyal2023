#####################NOTE#####################
#FateMap Revisions
#First started by YG on March 14, 2022
#Barcode Levenstein Distance Cutoff analysis
##############################################

library(tidyverse)
library(stringdist)
library(ggplot2)
library(dplyr)
library(gridExtra)

#macDirectory <- '/Users/yogeshgoyal/'
#macDirectory <- '/Users/ygy1258/'

input1Directory <- 'Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM01/'
input2Directory <- 'Dropbox (RajLab)/FateMap/paper/rawData/10X_cellRangerOuts/FM01/1_1uMPLX/filtered_feature_bc_matrix/'
input3Directory<- "Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM01/"
plotDirectory <- "Dropbox (RajLab)/FateMap/REVISIONS/analysisPlots/"

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

######################Barcode 50
barcode50 = as_tibble(read.table(paste0(input1Directory, 'stepFourStarcodeShavedReads50.txt'), stringsAsFactors=F, header = T)) %>% filter(sampleNum == "1")
data1file = as_tibble(read.table(paste0(input2Directory,"barcodes.tsv.gz"), stringsAsFactors=F)) %>% dplyr::rename(cellID = V1) ###original barcode file
data1file = as_tibble(substring(data1file$cellID, 1,nchar(data1file[1,1])-2)) %>% dplyr::rename(cellID = value) 

umapCoordinatesS1 = as_tibble(read.table(file = paste0(input3Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% filter(sampleNum == "S1")
barcode50CellID = as_tibble(unique(barcode50$cellID)) %>% rename(cellID = value) 
umapNoCutoff50 = inner_join(barcode50CellID, umapCoordinatesS1, by = "cellID")

uniqueCellIDDetectedS1 = inner_join(barcode50CellID, umapCoordinatesS1, by = "cellID")
UnDetectedS1 = umapCoordinatesS1 %>% filter(!cellID %in% uniqueCellIDDetectedS1$cellID)

initialPercentDetectedS1 = 100*length(unique(uniqueCellIDDetectedS1$cellID))/length(unique(umapCoordinatesS1$cellID))

plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray50", size = 1.5, shape = 16) +
  geom_point(data = uniqueCellIDDetectedS1, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3", size = 1.5, shape = 16, alpha = 1) +
  #geom_point(data = UnDetectedS1, aes(x = UMAP_1, y = UMAP_2), color = "gray94", size = 1.5, shape = 16, alpha = 0.5) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plotDirectory, 'FM01_initialPercentDetected.svg'), width = 6, height = 5.325)

umiCut = 4 # Minimum UMI cutoff for reliable analysis. Can be played around with between 3 to 5.
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.

upToLineageCounting = barcode50 %>% 
  filter(sampleNum == "1") %>%
  group_by(cellID, BC50StarcodeD8, sampleNum) %>%
  summarise(nUMI = length(sampleNum)) %>%
  filter(nUMI >= umiCut) %>% 
  group_by(cellID) %>%
  mutate(nLineages = length(cellID)) 

uniqueCellIDUpToLineageCounting = inner_join((as_tibble(unique(upToLineageCounting$cellID)) %>% rename(cellID = value)), umapCoordinatesS1 )
  
UMICutoffPercentDetectedS1 = 100*length(unique(uniqueCellIDUpToLineageCounting$cellID))/length(unique(uniqueCellIDDetectedS1$cellID))
upToLineageCounting$sampleNum = sub("^", "S", upToLineageCounting$sampleNum)

upToLineageCountingSummary = upToLineageCounting %>% select(-BC50StarcodeD8, -nUMI) %>% unique() %>% group_by(nLineages) %>% summarise(summaryLineages = length(nLineages))

write.table(upToLineageCountingSummary, file=paste0(plotDirectory,'summaryBarcodeCells.tsv'), col.names = TRUE, sep='\t')

plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray50", size = 1.5, shape = 16) +
  geom_point(data = uniqueCellIDUpToLineageCounting, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3", size = 1.5, shape = 16, alpha = 1) +
  #geom_point(data = UnDetectedS1, aes(x = UMAP_1, y = UMAP_2), color = "gray94", size = 1.5, shape = 16, alpha = 0.5) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plotDirectory, 'FM01_uniqueCellIDUpToLineageCounting.svg'), width = 6, height = 5.325)

#####
###taking only single cellid-->barcode mappings
linCountTooverlaps = upToLineageCounting %>%
  ungroup() %>%
  filter(nLineages <= linCut) %>%
  unique()


uniqueCellIDLinCountTooverlaps = inner_join((as_tibble(unique(linCountTooverlaps$cellID)) %>% rename(cellID = value)), umapCoordinatesS1, by = "cellID")
finalPercentDetectedPostFiltering50 = 100*length(unique(uniqueCellIDLinCountTooverlaps$cellID))/length(unique(uniqueCellIDUpToLineageCounting$cellID))

plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray50", size = 1.5, shape = 16) +
  geom_point(data = uniqueCellIDLinCountTooverlaps, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3", size = 1.5, shape = 16, alpha = 1) +
  #geom_point(data = UnDetectedS1, aes(x = UMAP_1, y = UMAP_2), color = "gray94", size = 1.5, shape = 16, alpha = 0.5) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plotDirectory, 'uniqueCellIDLinCountTooverlaps.svg'), width = 6, height = 5.325)


####doublet calculations
# for 0.23 MOI = (0.20546-0.18274)/0.20546 = ~11.05%
# Doublets from 10X experiment = ~10%
#Total unexplained = 29-11-10 = ~8%

####Clone Size Distribution Calculations

jointUMAP = inner_join(linCountTooverlaps, umapCoordinatesS1, by = c("cellID", "sampleNum")) %>% select(-nLineages)
colonyDistributionS150 = jointUMAP %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) 
colonyDistributionS150NonSinglet =colonyDistributionS150 %>% filter(nColony >1)

cloneDistribution = ggplot(colonyDistributionS150, aes(nColony)) +
  geom_bar() +
  geom_rug() +
  geom_vline(xintercept= mean(colonyDistributionS150$nColony), color="blue", size=0.5) + #2.27
  geom_vline(xintercept= median(colonyDistributionS150$nColony), color="red", size=0.5) +
  ylim(0,2000) +
  theme_classic()
ggsave(cloneDistribution, file = paste0(plotDirectory, 'FM01_colonyDistributionS150.svg'), width = 6, height = 5.325)

cloneDistributionNonSinglet = ggplot(colonyDistributionS150NonSinglet, aes(nColony)) +
  geom_bar() +
  geom_rug() +
  geom_vline(xintercept= mean(colonyDistributionS150NonSinglet$nColony), color="blue", size=0.5) + #12.45
  ylim(0,200) +
  theme_classic()
ggsave(cloneDistributionNonSinglet, file = paste0(plotDirectory, 'FM01_colonyDistributionS150NonSinglet.svg'), width = 6, height = 5.325)

######################Barcode 40
barcode40 = as_tibble(read.table(paste0(input1Directory, 'stepFourStarcodeShavedReads40.txt'), stringsAsFactors=F, header = T))
umiCut = 4 # Minimum UMI cutoff for reliable analysis. Can be played around with between 3 to 5.
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.  

uniqueCellIDDetectedS1 = barcode40 %>% filter(sampleNum == 1)

initialPercentDetected = 100*length(unique(uniqueCellIDDetectedS1$cellID))/length(unique(data1file$cellID))

upToLineageCounting = barcode40 %>% 
  filter(sampleNum == "1") %>%
  group_by(cellID, BC40StarcodeD8, sampleNum) %>%
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

finalPercentDetectedPostFiltering40 = 100*length(unique(linCountTooverlaps$cellID))/length(unique(uniqueCellIDDetectedS1$cellID))

####Clone Size Distribution Calculations

jointUMAP = inner_join(linCountTooverlaps, umapCoordinatesS1, by = c("cellID", "sampleNum")) %>% select(-nLineages)
colonyDistributionS140 = jointUMAP %>% group_by(BC40StarcodeD8) %>% summarise(nColony = length(BC40StarcodeD8)) 
colonyDistributionS140NonSinglet =colonyDistributionS140 %>% filter(nColony >1)

cloneDistribution = ggplot(colonyDistributionS140, aes(nColony)) +
  geom_bar() +
  geom_rug() +
  geom_vline(xintercept= mean(colonyDistributionS140$nColony), color="blue", size=0.5) + #2.27
  geom_vline(xintercept= median(colonyDistributionS140$nColony), color="red", size=0.5) +
  ylim(0,2000) +
  theme_classic()
ggsave(cloneDistribution, file = paste0(plotDirectory, 'FM01_colonyDistributionS140.svg'), width = 6, height = 5.325)

cloneDistributionNonSinglet = ggplot(colonyDistributionS140NonSinglet, aes(nColony)) +
  geom_bar() +
  geom_rug() +
  geom_vline(xintercept= mean(colonyDistributionS140NonSinglet$nColony), color="blue", size=0.5) + #12.45
  ylim(0,200) +
  theme_classic()
ggsave(cloneDistributionNonSinglet, file = paste0(plotDirectory, 'FM01_colonyDistributionS140NonSinglet.svg'), width = 6, height = 5.325)

######################Barcode 40
barcode30 = as_tibble(read.table(paste0(input1Directory, 'stepFourStarcodeShavedReads30.txt'), stringsAsFactors=F, header = T))
umiCut = 4 # Minimum UMI cutoff for reliable analysis. Can be played around with between 3 to 5.
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.  

uniqueCellIDDetectedS1 = barcode30 %>% filter(sampleNum == 1)

initialPercentDetected = 100*length(unique(uniqueCellIDDetectedS1$cellID))/length(unique(data1file$cellID))

upToLineageCounting = barcode30 %>% 
  filter(sampleNum == "1") %>%
  group_by(cellID, BC30StarcodeD8, sampleNum) %>%
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

finalPercentDetectedPostFiltering30 = 100*length(unique(linCountTooverlaps$cellID))/length(unique(uniqueCellIDDetectedS1$cellID))


jointUMAP = inner_join(linCountTooverlaps, umapCoordinatesS1, by = c("cellID", "sampleNum")) %>% select(-nLineages)
colonyDistributionS130 = jointUMAP %>% group_by(BC30StarcodeD8) %>% summarise(nColony = length(BC30StarcodeD8)) 
colonyDistributionS130NonSinglet =colonyDistributionS130 %>% filter(nColony >1)

cloneDistribution = ggplot(colonyDistributionS130, aes(nColony)) +
  geom_bar() +
  geom_rug() +
  geom_vline(xintercept= mean(colonyDistributionS130$nColony), color="blue", size=0.5) + #2.27
  geom_vline(xintercept= median(colonyDistributionS130$nColony), color="red", size=0.5) +
  ylim(0,2000) +
  theme_classic()
ggsave(cloneDistribution, file = paste0(plotDirectory, 'FM01_colonyDistributionS130.svg'), width = 6, height = 5.325)

cloneDistributionNonSinglet = ggplot(colonyDistributionS130NonSinglet, aes(nColony)) +
  geom_bar() +
  geom_rug() +
  geom_vline(xintercept= mean(colonyDistributionS130NonSinglet$nColony), color="blue", size=0.5) + #12.45
  ylim(0,200) +
  theme_classic()
ggsave(cloneDistributionNonSinglet, file = paste0(plotDirectory, 'FM01_colonyDistributionS130NonSinglet.svg'), width = 6, height = 5.325)

########
## Robustness to UMI cutoff for LV50d8

umiCutRange = c(1:15)
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.
uniqueCellIDUpToLineageCountingTally = tibble(UMICutoff = numeric(),
                                              totalCellsUMI = numeric(),
                                              totalCellsLinCut = numeric(),
                                              colonyMean = numeric())
for (umiCut in umiCutRange) {

  upToLineageCounting = barcode50 %>% 
    filter(sampleNum == "1") %>%
    group_by(cellID, BC50StarcodeD8, sampleNum) %>%
    summarise(nUMI = length(sampleNum)) %>%
    filter(nUMI >= umiCut) %>% 
    group_by(cellID) %>%
    mutate(nLineages = length(cellID)) 
  
  uniqueCellIDUpToLineageCounting = inner_join((as_tibble(unique(upToLineageCounting$cellID)) %>% rename(cellID = value)), umapCoordinatesS1, by = "cellID")
  uniqueCellIDUpToLineageCountingTally[umiCut,1] = umiCut - 1
  uniqueCellIDUpToLineageCountingTally[umiCut,2] = length(unique(uniqueCellIDUpToLineageCounting$cellID))
  #####
  ###taking only single cellid-->barcode mappings
  linCountTooverlaps = upToLineageCounting %>%
    ungroup() %>%
    filter(nLineages <= linCut) %>%
    unique()

  linCountTooverlaps$sampleNum = sub("^", "S", linCountTooverlaps$sampleNum) 

  uniqueCellIDLinCountTooverlaps = inner_join((as_tibble(unique(linCountTooverlaps$cellID)) %>% rename(cellID = value)), umapCoordinatesS1, by = "cellID")
  uniqueCellIDUpToLineageCountingTally[umiCut,3] = length(unique(uniqueCellIDLinCountTooverlaps$cellID))
  
  jointUMAP = inner_join(linCountTooverlaps, umapCoordinatesS1, by = c("cellID", "sampleNum")) %>% select(-nLineages)
  colonyDistributionS150 = jointUMAP %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8))
  uniqueCellIDUpToLineageCountingTally[umiCut,4] = mean(colonyDistributionS150$nColony)
  
}

totalCells = ggplot(uniqueCellIDUpToLineageCountingTally, aes(UMICutoff,totalCellsUMI)) +
  geom_col()+
  scale_x_discrete(limits=c(0:14)) +
  theme_classic()
ggsave(totalCells, file = paste0(plotDirectory, 'FM01_totalCellsUMIOnly.svg'), width = 3, height = 2.7)


totalCells = ggplot(uniqueCellIDUpToLineageCountingTally, aes(UMICutoff,totalCellsLinCut)) +
  geom_col()+
  scale_x_discrete(limits=c(0:14)) +
  theme_classic()
ggsave(totalCells, file = paste0(plotDirectory, 'FM01_totalCellsUMILin.svg'), width = 3, height = 2.7)

cloneSize = ggplot(uniqueCellIDUpToLineageCountingTally, aes(UMICutoff,colonyMean)) +
  geom_col()+
  scale_x_discrete(limits=c(0:14)) +
  theme_classic()
ggsave(cloneSize, file = paste0(plotDirectory, 'FM01_cloneSize.svg'), width = 3, height = 2.7)


