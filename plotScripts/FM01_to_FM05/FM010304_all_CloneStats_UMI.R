#####################NOTE#####################
#FateMap Revisions_ 
#First started by YG on March 21, 2022
###last edited by YG Feb 05, 2023
#Barcode Recovery Across Cell Types
##############################################

library(tidyverse)
library(stringdist)
library(ggplot2)
library(dplyr)
library(gridExtra)

#macDirectory <- '/Users/yogeshgoyal/'
#macDirectory <- '/Users/ygy1258/'



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
#####WM989, FM01
input1Directory <- 'Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM01/'
input2Directory<- "Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM01/"
plotDirectory <- "Dropbox (RajLab)/FateMap/REVISIONS/analysisPlots/workingFolder/"

#####TESTING for Maalavika analysis
#testDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM01/"
#testBarcode50 = as_tibble(read.table(paste0(testDirectory, 'barcodeCellID.tsv'), stringsAsFactors=F, header = T)) 
#nrow(unique(testBarcode50 %>% filter(sampleNum =="S1")))
#nrow(testBarcode50 %>% filter(sampleNum =="S2"))
##############################################################

barcode50 = as_tibble(read.table(paste0(input1Directory, 'stepFourStarcodeShavedReads50.txt'), stringsAsFactors=F, header = T)) %>% filter(sampleNum %in% c("1","2"))

umapCoordinatesS1S2 = as_tibble(read.table(file = paste0(input2Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% filter(sampleNum %in% c("S1","S2"))
barcode50CellID = as_tibble(unique(barcode50$cellID)) %>% rename(cellID = value) 
uniqueCellIDDetectedS1S2 = inner_join(barcode50CellID, umapCoordinatesS1S2, by = "cellID")

initialPercentDetectedS1 = 100*nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S1"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S1"))
initialPercentDetectedS2 = 100*nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S2"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S2"))

umiCut = 4 # Minimum UMI cutoff for reliable analysis. Can be played around with between 3 to 5.
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.

upToLineageCounting = barcode50 %>% 
  filter(sampleNum %in% c("1","2")) %>%
  group_by(cellID, BC50StarcodeD8, sampleNum) %>%
  summarise(nUMI = length(sampleNum)) %>%
  filter(nUMI >= umiCut) %>% 
  group_by(cellID) %>%
  mutate(nLineages = length(cellID)) 

uniqueCellIDUpToLineageCounting = inner_join((as_tibble(unique(upToLineageCounting$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2)

UMICutoffPercentDetectedS1 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S1"))/nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S1"))
UMICutoffPercentDetectedS2 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S2"))/nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S2"))

upToLineageCounting$sampleNum = sub("^", "S", upToLineageCounting$sampleNum)

###taking only single cellid-->barcode mappings
linCountTooverlaps = upToLineageCounting %>%
  ungroup() %>%
  filter(nLineages <= linCut) %>%
  unique()

uniqueCellIDLinCountTooverlaps = inner_join((as_tibble(unique(linCountTooverlaps$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2, by = "cellID")

nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S1"))
nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S2"))

##S1 = 4672; S2 = 4452; Date added: Feb 05, 2023

finalPercentDetectedPostFilteringS1 = 100*nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S1"))/nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S1"))
finalPercentDetectedPostFilteringS2 = 100*nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S2"))/nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S2"))


#####WM989, FM02 S1,S2
input1Directory <- 'Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM02/'
input2Directory<- "Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM02/"

barcode50 = as_tibble(read.table(paste0(input2Directory, 'stepThreeStarcodeShavedReads.txt'), stringsAsFactors=F, header = F))
barcode50 = barcode50  %>% 
  dplyr::rename(cellID = V1,
                UMI = V2,
                BC50StarcodeD8 = V4,
                sampleNum = V8) %>%
  select(-c(V3,V5,V6,V7)) %>% unique() %>% select(-UMI) %>% filter(sampleNum %in% c("1","2", "3", "4"))

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

nrow(linCountTooverlaps %>% filter(sampleNum =="S1"))
nrow(linCountTooverlaps %>% filter(sampleNum =="S2"))
nrow(linCountTooverlaps %>% filter(sampleNum =="S3"))
nrow(linCountTooverlaps %>% filter(sampleNum =="S4"))

##S1 = 5823; S2 = 6232; S3:5193 S4: 4447

umapCoordinatesS1S2 = as_tibble(read.table(file = paste0(input1Directory, "umapCoordinates_s1s2Scanorama_50pcs_filter_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% filter(sampleNum %in% c("S1","S2", "S3"))
umapCoordinatesS1S2 = as_tibble(read.table(file = paste0(input1Directory, "umapCoordinates_s1s3s4Scanorama_50pcs_filter_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% filter(sampleNum %in% c("S1","S3", "S4"))

barcode50CellID = as_tibble(unique(barcode50$cellID)) %>% rename(cellID = value) 

uniqueCellIDDetectedS1S2 = inner_join(barcode50CellID, umapCoordinatesS1S2, by = "cellID")

initialPercentDetectedS1 = 100*nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S1"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S1"))
initialPercentDetectedS2 = 100*nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S2"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S2"))

uniqueCellIDUpToLineageCounting = inner_join((as_tibble(unique(upToLineageCounting$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2)

UMICutoffPercentDetectedS1 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S1"))/nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S1"))
UMICutoffPercentDetectedS2 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S2"))/nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S2"))

upToLineageCounting$sampleNum = sub("^", "S", upToLineageCounting$sampleNum) 
#####

uniqueCellIDLinCountTooverlaps = inner_join((as_tibble(unique(linCountTooverlaps$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2, by = "cellID")

nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S1"))
nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S3"))
nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S4"))

##S1 = 5826; S2 = 6223; S3: 5200 ; S4: 4447  Date added: Feb 05, 2023

#####WM989, FM03 S1,S2

input1Directory <- 'Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM03/'
input2Directory<- "Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM03/"

barcode50 = as_tibble(read.table(paste0(input2Directory, 'stepThreeStarcodeShavedReads.txt'), stringsAsFactors=F, header = F))
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

nrow(linCountTooverlaps %>% filter(sampleNum =="S1"))
nrow(linCountTooverlaps %>% filter(sampleNum =="S2"))

##S1 = 4351; S2 = 3854

umapCoordinatesS1S2 = as_tibble(read.table(file = paste0(input1Directory, "umapCoordinates_s1s2Scanorama_50pcs_filter_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% filter(sampleNum %in% c("S1","S2"))

barcode50CellID = as_tibble(unique(barcode50$cellID)) %>% rename(cellID = value) 

uniqueCellIDDetectedS1S2 = inner_join(barcode50CellID, umapCoordinatesS1S2, by = "cellID")

initialPercentDetectedS1 = 100*nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S1"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S1"))
initialPercentDetectedS2 = 100*nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S2"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S2"))

uniqueCellIDUpToLineageCounting = inner_join((as_tibble(unique(upToLineageCounting$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2)

UMICutoffPercentDetectedS1 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S1"))/nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S1"))
UMICutoffPercentDetectedS2 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S2"))/nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S2"))

upToLineageCounting$sampleNum = sub("^", "S", upToLineageCounting$sampleNum) 
#####

uniqueCellIDLinCountTooverlaps = inner_join((as_tibble(unique(linCountTooverlaps$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2, by = "cellID")

nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S1"))
nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S2"))

##S1 = 4354; S2 = 3842; Date added: Feb 05, 2023

#####WM983B, FM05
input1Directory <- 'Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM05/'
input2Directory<- "Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM05/"

barcode50 = as_tibble(read.table(paste0(input2Directory, 'stepThreeStarcodeShavedReads.txt'), stringsAsFactors=F, header = F)) 
barcode50 = barcode50  %>% 
  dplyr::rename(cellID = V1,
                UMI = V2,
                BC50StarcodeD8 = V4,
                sampleNum = V8) %>%
  select(-c(V3,V5,V6,V7)) %>% unique() %>% select(-UMI) %>% filter(sampleNum %in% c("1","2","3","4"))

umapCoordinatesS1S2 = as_tibble(read.table(file = paste0(input1Directory, "umapCoordinatesSCT_S1S2_FM05_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% filter(sampleNum %in% c("S1","S2"))
barcode50CellID = as_tibble(unique(barcode50$cellID)) %>% rename(cellID = value) 

uniqueCellIDDetectedS1S2 = inner_join(barcode50CellID, umapCoordinatesS1S2, by = "cellID")

initialPercentDetectedS1 = 100*nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S1"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S1"))
initialPercentDetectedS2 = 100*nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S2"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S2"))

umiCut = 3 # Minimum UMI cutoff for reliable analysis. For the past runs I used between 3-5 for ~500k reads/sample.
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.  

upToLineageCounting = barcode50 %>% 
  group_by(cellID, BC50StarcodeD8, sampleNum) %>%
  summarise(nUMI = length(sampleNum)) %>%
  filter(nUMI >= umiCut) %>% 
  group_by(cellID) %>%
  mutate(nLineages = length(cellID)) 

uniqueCellIDUpToLineageCounting = inner_join((as_tibble(unique(upToLineageCounting$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2)

UMICutoffPercentDetectedS1 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S1"))/nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S1"))
UMICutoffPercentDetectedS2 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S2"))/nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S2"))

upToLineageCounting$sampleNum = sub("^", "S", upToLineageCounting$sampleNum) 
#####
###taking only single cellid-->barcode mappings
linCountTooverlaps = upToLineageCounting %>%
  ungroup() %>%
  filter(nLineages <= linCut) %>%
  unique()

nrow(linCountTooverlaps %>% filter(sampleNum =="S1"))
nrow(linCountTooverlaps %>% filter(sampleNum =="S2"))

##S1 = 5553; S2 = 6475

uniqueCellIDLinCountTooverlaps = inner_join((as_tibble(unique(linCountTooverlaps$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2, by = "cellID")

nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S1"))
nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S2"))

##S1 = 5553; S2 = 6475; Date added: Feb 05, 2023

finalPercentDetectedPostFilteringS1 = 100*nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S1"))/nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S1"))
finalPercentDetectedPostFilteringS2 = 100*nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S2"))/nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S2"))


#####MDA, FM04
input1Directory <- 'Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM04/'
input2Directory<- "Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM04/"

barcode50 = as_tibble(read.table(paste0(input1Directory, 'stepThreeStarcodeShavedReads.txt'), stringsAsFactors=F, header = F)) 
barcode50 = barcode50  %>% 
  dplyr::rename(cellID = V1,
                UMI = V2,
                BC50StarcodeD8 = V4,
                sampleNum = V8) %>%
  select(-c(V3,V5,V6,V7)) %>% unique() %>% select(-UMI) %>% filter(sampleNum %in% c("1","2"))

umapCoordinatesS1S2 = as_tibble(read.table(file = paste0(input2Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% filter(sampleNum %in% c("S1","S2"))
barcode50CellID = as_tibble(unique(barcode50$cellID)) %>% rename(cellID = value) 

uniqueCellIDDetectedS1S2 = inner_join(barcode50CellID, umapCoordinatesS1S2, by = "cellID")

initialPercentDetected = 100*nrow(uniqueCellIDDetectedS1S2 )/nrow(umapCoordinatesS1S2)


initialPercentDetectedS1 = 100*nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S1"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S1"))
initialPercentDetectedS2 = 100*nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S2"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S2"))

umiCut = 4 # Minimum UMI cutoff for reliable analysis. For the past runs I used between 3-5 for ~500k reads/sample.
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.  

upToLineageCounting = barcode50 %>% 
  group_by(cellID, BC50StarcodeD8, sampleNum) %>%
  summarise(nUMI = length(sampleNum)) %>%
  filter(nUMI >= umiCut) %>% 
  group_by(cellID) %>%
  mutate(nLineages = length(cellID)) 

uniqueCellIDUpToLineageCounting = inner_join((as_tibble(unique(upToLineageCounting$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2)

UMICutoffPercentDetectedS1 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S1"))/nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S1"))
UMICutoffPercentDetectedS2 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S2"))/nrow(uniqueCellIDDetectedS1S2 %>% filter(sampleNum =="S2"))

plot <- ggplot() +
  geom_point(data = umapCoordinatesS1S2, aes(x = UMAP_1, y = UMAP_2), color = "gray55", size = 1.5, shape = 16) +
  geom_point(data = uniqueCellIDUpToLineageCounting, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3", size = 1.5, shape = 16, alpha = 1) +
  #geom_point(data = UnDetectedS1, aes(x = UMAP_1, y = UMAP_2), color = "gray94", size = 1.5, shape = 16, alpha = 0.5) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plotDirectory, 'FM04_InitialPercentDetected.svg'), width = 6, height = 5)

upToLineageCounting$sampleNum = sub("^", "S", upToLineageCounting$sampleNum) 
#####
###taking only single cellid-->barcode mappings
linCountTooverlaps = upToLineageCounting %>%
  ungroup() %>%
  filter(nLineages <= linCut) %>%
  unique()

nrow(linCountTooverlaps %>% filter(sampleNum =="S1"))
nrow(linCountTooverlaps %>% filter(sampleNum =="S2"))

##S1 = 3174; S2 = 4665

uniqueCellIDLinCountTooverlaps = inner_join((as_tibble(unique(linCountTooverlaps$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2, by = "cellID")

nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S1"))
nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S2"))

##S1 = 5460; S2 = 6367; Date added: Feb 05, 2023

finalPercentDetectedPostFilteringS1 = 100*nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S1"))/nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S1"))
finalPercentDetectedPostFilteringS2 = 100*nrow(uniqueCellIDLinCountTooverlaps %>% filter(sampleNum =="S2"))/nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S2"))

plot <- ggplot() +
  geom_point(data = umapCoordinatesS1S2, aes(x = UMAP_1, y = UMAP_2), color = "gray55", size = 1.5, shape = 16) +
  geom_point(data = uniqueCellIDLinCountTooverlaps, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3", size = 1.5, shape = 16, alpha = 1) +
  #geom_point(data = UnDetectedS1, aes(x = UMAP_1, y = UMAP_2), color = "gray94", size = 1.5, shape = 16, alpha = 0.5) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plotDirectory, 'FM04_usedPercentDetected.svg'), width = 6, height = 5)

#####
####FM06, WM989 Naive
input2Directory<- 'Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2021_FM06/extracted/FM06_barcode/'
barcodesCellIDs = as_tibble(read.table(paste0(input2Directory, 'barcodeCellID.tsv'), stringsAsFactors=F, header = T))
nrow(barcodesCellIDs %>% filter(sampleNum =="S1"))
nrow(barcodesCellIDs %>% filter(sampleNum =="S2"))
                  
# barcode50 = as_tibble(read.table(paste0(input2Directory, 'stepThreeStarcodeShavedReads.txt'), stringsAsFactors=F, header = F))
# barcode50 = barcode50  %>%  
#   dplyr::rename(cellID = V1,
#                 UMI = V2,
#                 BC50StarcodeD8 = V4,
#                 sampleNum = V8) %>%
#   select(-c(V3,V5,V6,V7)) %>% unique() %>% select(-UMI) %>% filter(sampleNum %in% c("1","2"))
# 
# umiCut = 4 # Minimum UMI cutoff for reliable analysis. For the past runs I used between 3-5 for ~500k reads/sample. Here we have ~1.5mil/sample
# linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.  
# 
# upToLineageCounting = barcode50 %>% 
#   group_by(cellID, BC50StarcodeD8, sampleNum) %>%
#   summarise(nUMI = length(sampleNum)) %>%
#   filter(nUMI >= umiCut) %>% 
#   group_by(cellID) %>%
#   mutate(nLineages = length(cellID)) 
# 
# upToLineageCounting$sampleNum = sub("^", "S", upToLineageCounting$sampleNum) 
# 
# ###taking only single cellid-->barcode mappings
# linCountTooverlaps = upToLineageCounting %>%
#   ungroup() %>%
#   filter(nLineages <= linCut) %>%
#   unique()

####FM07, WM989 Cont vs Discont

input2Directory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM07/jointAnalysis/seurat/"
barcodesCellIDs <- read.table(paste0(input2Directory,'barcodeCellID2.tsv'), stringsAsFactors=F, header = T, sep="\t")
nrow(barcodesCellIDs %>% filter(sampleNum =="S1"))
nrow(barcodesCellIDs %>% filter(sampleNum =="S2"))

####FM08, Melanocytes

input2Directory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM09/jointAnalysis/seurat/"
barcodesCellIDs <- read.table(paste0(input2Directory,'barcodeCellID.tsv'), stringsAsFactors=F, header = T, sep="\t")
nrow(barcodesCellIDs %>% filter(sampleNum =="S1"))
nrow(barcodesCellIDs %>% filter(sampleNum =="S2"))

####FM09, NRAS mutant WM3451

input2Directory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM09/jointAnalysis/seurat/"

barcodesCellIDs <- read.table(paste0(input2Directory,'barcodeCellID.tsv'), stringsAsFactors=F, header = T, sep="\t")
nrow(barcodesCellIDs %>% filter(sampleNum =="S1"))
nrow(barcodesCellIDs %>% filter(sampleNum =="S2"))

####FM09, NRAS mutant WM

input2Directory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM10/jointAnalysis/seurat/"

barcodesCellIDs <- read.table(paste0(input2Directory,'barcodeCellID.tsv'), stringsAsFactors=F, header = T, sep="\t")
nrow(barcodesCellIDs %>% filter(sampleNum =="S1"))
nrow(barcodesCellIDs %>% filter(sampleNum =="S2"))

###############################################33
####FOR DIFFERENT CUTOFFS
######################
input1Directory <- 'Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM01/'
input2Directory<- "Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM01/"
plotDirectory <- "Dropbox (RajLab)/FateMap/REVISIONS/analysisPlots/"

barcode50 = as_tibble(read.table(paste0(input1Directory, 'stepFourStarcodeShavedReads50.txt'), stringsAsFactors=F, header = T)) %>% filter(sampleNum %in% c("1","2"))
umapCoordinatesS1S2 = as_tibble(read.table(file = paste0(input2Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% filter(sampleNum %in% c("S1","S2"))

umiCutRange = c(1:10)
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.
uniqueCellIDUpToLineageCountingTally = tibble(UMICutoff = numeric(),
                                              percentS1 = numeric(),
                                              percentS2 = numeric(),
                                              cellType = character())
for (umiCut in umiCutRange) {
  
  upToLineageCounting = barcode50 %>% 
    filter(sampleNum %in% c("1","2")) %>%
    group_by(cellID, BC50StarcodeD8, sampleNum) %>%
    summarise(nUMI = length(sampleNum)) %>%
    filter(nUMI >= umiCut) %>% 
    group_by(cellID) %>%
    mutate(nLineages = length(cellID)) 

  uniqueCellIDUpToLineageCounting = inner_join((as_tibble(unique(upToLineageCounting$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2)
  
  initialPercentDetectedS1 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S1"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S1"))
  initialPercentDetectedS2 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S2"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S2"))
  
  uniqueCellIDUpToLineageCountingTally[umiCut,1] = umiCut - 1
  uniqueCellIDUpToLineageCountingTally[umiCut,2] = initialPercentDetectedS1
  uniqueCellIDUpToLineageCountingTally[umiCut,3] = initialPercentDetectedS2
  uniqueCellIDUpToLineageCountingTally[umiCut,4] = "wm989"
}

######################3
#####WM983B, FM05
input1Directory <- 'Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM05/'
input2Directory<- "Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM05/"

barcode50 = as_tibble(read.table(paste0(input2Directory, 'stepThreeStarcodeShavedReads.txt'), stringsAsFactors=F, header = F)) 
barcode50 = barcode50  %>% 
  dplyr::rename(cellID = V1,
                UMI = V2,
                BC50StarcodeD8 = V4,
                sampleNum = V8) %>%
  select(-c(V3,V5,V6,V7)) %>% unique() %>% select(-UMI) %>% filter(sampleNum %in% c("1","2"))

umapCoordinatesS1S2 = as_tibble(read.table(file = paste0(input1Directory, "umapCoordinatesSCT_S1S2_FM05_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% filter(sampleNum %in% c("S1","S2"))

umiCutRange = c(1:10)
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.

for (umiCut in umiCutRange) {
  
  upToLineageCounting = barcode50 %>% 
    filter(sampleNum %in% c("1","2")) %>%
    group_by(cellID, BC50StarcodeD8, sampleNum) %>%
    summarise(nUMI = length(sampleNum)) %>%
    filter(nUMI >= umiCut) %>% 
    group_by(cellID) %>%
    mutate(nLineages = length(cellID)) 
  
  uniqueCellIDUpToLineageCounting = inner_join((as_tibble(unique(upToLineageCounting$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2)
  
  initialPercentDetectedS1 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S1"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S1"))
  initialPercentDetectedS2 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S2"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S2"))
  
  uniqueCellIDUpToLineageCountingTally[umiCut+10,1] = umiCut - 1
  uniqueCellIDUpToLineageCountingTally[umiCut+10,2] = initialPercentDetectedS1
  uniqueCellIDUpToLineageCountingTally[umiCut+10,3] = initialPercentDetectedS2
  uniqueCellIDUpToLineageCountingTally[umiCut+10,4] = "wm983b"
}

#####
#####MDA, FM04
input1Directory <- 'Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM04/'
input2Directory<- "Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM04/"

barcode50 = as_tibble(read.table(paste0(input1Directory, 'stepThreeStarcodeShavedReads.txt'), stringsAsFactors=F, header = F)) 
barcode50 = barcode50  %>% 
  dplyr::rename(cellID = V1,
                UMI = V2,
                BC50StarcodeD8 = V4,
                sampleNum = V8) %>%
  select(-c(V3,V5,V6,V7)) %>% unique() %>% select(-UMI) %>% filter(sampleNum %in% c("1","2"))

umapCoordinatesS1S2 = as_tibble(read.table(file = paste0(input2Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% filter(sampleNum %in% c("S1","S2"))

umiCutRange = c(1:10)
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.

for (umiCut in umiCutRange) {
  
  upToLineageCounting = barcode50 %>% 
    filter(sampleNum %in% c("1","2")) %>%
    group_by(cellID, BC50StarcodeD8, sampleNum) %>%
    summarise(nUMI = length(sampleNum)) %>%
    filter(nUMI >= umiCut) %>% 
    group_by(cellID) %>%
    mutate(nLineages = length(cellID)) 
  
  uniqueCellIDUpToLineageCounting = inner_join((as_tibble(unique(upToLineageCounting$cellID)) %>% rename(cellID = value)), umapCoordinatesS1S2)
  
  initialPercentDetectedS1 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S1"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S1"))
  initialPercentDetectedS2 = 100*nrow(uniqueCellIDUpToLineageCounting %>% filter(sampleNum =="S2"))/nrow(umapCoordinatesS1S2 %>% filter(sampleNum =="S2"))
  
  uniqueCellIDUpToLineageCountingTally[umiCut+20,1] = umiCut - 1
  uniqueCellIDUpToLineageCountingTally[umiCut+20,2] = initialPercentDetectedS1
  uniqueCellIDUpToLineageCountingTally[umiCut+20,3] = initialPercentDetectedS2
  uniqueCellIDUpToLineageCountingTally[umiCut+20,4] = "MDA"
}

uniqueCellIDUpToLineageCountingTally = uniqueCellIDUpToLineageCountingTally %>% mutate(meanPercentRecovered = (percentS1+percentS2)/2) %>% mutate(UMICutoff = UMICutoff + 1)

ggplot(uniqueCellIDUpToLineageCountingTally, aes(UMICutoff, meanPercentRecovered)) +
  geom_point()+
  geom_line()+
  facet_grid(cols = vars(cellType)) +
  theme_classic()


