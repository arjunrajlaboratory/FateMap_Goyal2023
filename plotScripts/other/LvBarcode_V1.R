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

macDirectory <- '/Users/yogeshgoyal/'
macDirectory <- '/Users/ygy1258/'

input1Directory <- 'Dropbox (RajLab)/FateMap/paper/rawData/10X_cellRangerOuts/FM01/1_1uMPLX/filtered_feature_bc_matrix/'
input1aDirectory <- 'Dropbox (RajLab)/FateMap/paper/rawData/10X_cellRangerOuts/FM01/2_1uMPLX/filtered_feature_bc_matrix/'
input2Directory <- 'Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM01/1_1uMPLX/'
input2aDirectory <- 'Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM01/2_1uMPLX/'
input3Directory <- 'Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM01/'

plotDirectory <- "Dropbox (RajLab)/FateMap/REVISIONS/analysisPlots/"

data1file = as_tibble(read.table(paste0(macDirectory, input1Directory,"barcodes.tsv.gz"), stringsAsFactors=F)) %>% dplyr::rename(cellID = V1) 
data1file = as_tibble(substring(data1file$cellID, 1,nchar(data1file[1,1])-2)) %>% dplyr::rename(cellID = value) 
data2file = as_tibble(read.table(paste0(macDirectory, input2Directory,"uniqueShavedReads.txt"), stringsAsFactors=F)) %>% dplyr::rename(cellID = V1,UMI = V2, BC = V3) %>%
  mutate(BC50 = substring(BC,1,50),
         BC40 = substring(BC,1,40),
         BC30a = substring(BC,1,30),
         BC30b = substring(BC,1,30))
cellIDUMIBarcodes = inner_join(data1file, data2file, by = "cellID")
Barcodes = unique(cellIDUMIBarcodes$BC50)
cellIDs = unique(cellIDUMIBarcodes$cellID)

set.seed(2059)
subsample1 = sample(Barcodes,1000)
subsample2 = sample(Barcodes,1000)
subsample3 = sample(Barcodes,1000)
BarcodesLv1 = as.integer(stringdistmatrix(subsample1, method = "lv"))
BarcodesLv2 = as.integer(stringdistmatrix(subsample2, method = "lv"))
BarcodesLv3 = as.integer(stringdistmatrix(subsample3, method = "lv"))
lBarcodesLv = length(BarcodesLv1)

BarcodesLv = tibble(
  lvdist = c(BarcodesLv1, BarcodesLv2, BarcodesLv3),
  subsamNum = c(rep("subsamping1", lBarcodesLv), rep("subsamping2", lBarcodesLv), rep("subsamping3", lBarcodesLv)))

BarcodesLvHist <- BarcodesLv %>% group_by(subsamNum, lvdist) %>% summarise(length(lvdist)) %>%
  group_by(subsamNum) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

meanLvHisr = BarcodesLvHist %>% group_by(subsamNum) %>% summarise(mean = mean(lvdist)) #20, 19.5,19
ErrorMeanLvHisrS1 = meanLvHisr %>% mutate(errorPercent = (100*(25-mean)/25), condition = "S1Raw")

FM01_RawBarcodesLvHistPlot <- ggplot(BarcodesLvHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  facet_wrap(facets = vars(subsamNum)) +
  theme_classic() + 
  xlim(0,50) +
  geom_vline(aes(xintercept= mean),meanLvHisr, color="blue", size=0.5) +
  geom_vline(xintercept= 25, color="red", size=0.5)

FM01_RawBarcodesLvHistPlotMagnified <- ggplot(BarcodesLvHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  facet_wrap(facets = vars(subsamNum)) +
  theme_classic() + 
  xlim(0,30) +
  ylim(0,0.03)

set.seed(2059)
cellIDsLv = tibble(lvdist = as.integer(stringdistmatrix(cellIDs, method = "lv")))
cellIDsHist <- cellIDsLv  %>% group_by(lvdist)%>% summarise(length(lvdist)) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

FM01_cellIDsHistPlot <- ggplot(cellIDsHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  theme_classic()

ggsave(FM01_RawBarcodesLvHistPlot, file = paste0(plotDirectory, 'FM01_RawBarcodesLvHistPlotS1.svg'), width = 10, height = 4)
ggsave(FM01_cellIDsHistPlot, file = paste0(plotDirectory, 'FM01_cellIDsHistPlotS1.svg'), width = 4, height = 3)
ggsave(FM01_RawBarcodesLvHistPlotMagnified, file = paste0(plotDirectory, 'FM01_RawBarcodesLvHistPlotMagnifiedS1.svg'), width = 10, height = 4)

#####Sample 2 Raw
data1file = as_tibble(read.table(paste0(macDirectory, input1aDirectory,"barcodes.tsv.gz"), stringsAsFactors=F)) %>% dplyr::rename(cellID = V1) 
data1file = as_tibble(substring(data1file$cellID, 1,nchar(data1file[1,1])-2)) %>% dplyr::rename(cellID = value) 
data2file = as_tibble(read.table(paste0(macDirectory, input2aDirectory,"uniqueShavedReads.txt"), stringsAsFactors=F)) %>% dplyr::rename(cellID = V1,UMI = V2, BC = V3) %>%
  mutate(BC50 = substring(BC,1,50),
         BC40 = substring(BC,1,40),
         BC30a = substring(BC,1,30),
         BC30b = substring(BC,1,30))
cellIDUMIBarcodes = inner_join(data1file, data2file, by = "cellID")
Barcodes = unique(cellIDUMIBarcodes$BC50)
cellIDs = unique(cellIDUMIBarcodes$cellID)

set.seed(2059)
subsample1 = sample(Barcodes,1000)
subsample2 = sample(Barcodes,1000)
subsample3 = sample(Barcodes,1000)
BarcodesLv1 = as.integer(stringdistmatrix(subsample1, method = "lv"))
BarcodesLv2 = as.integer(stringdistmatrix(subsample2, method = "lv"))
BarcodesLv3 = as.integer(stringdistmatrix(subsample3, method = "lv"))
lBarcodesLv = length(BarcodesLv1)

BarcodesLv = tibble(
  lvdist = c(BarcodesLv1, BarcodesLv2, BarcodesLv3),
  subsamNum = c(rep("subsamping1", lBarcodesLv), rep("subsamping2", lBarcodesLv), rep("subsamping3", lBarcodesLv)))

BarcodesLvHist <- BarcodesLv %>% group_by(subsamNum, lvdist) %>% summarise(length(lvdist)) %>%
  group_by(subsamNum) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

meanLvHisr = BarcodesLvHist %>% group_by(subsamNum) %>% summarise(mean = mean(lvdist)) #20, 19.5,19
ErrorMeanLvHisrS2 = meanLvHisr %>% mutate(errorPercent = (100*(25-mean)/25), condition = "S2Raw")

############Once We ran Starcode

data3file = as_tibble(read.table(paste0(macDirectory, input3Directory, 'stepFourStarcodeShavedReads50.txt'), stringsAsFactors=F, header = T))

S1 = data3file %>% filter(sampleNum == 1)
S2 = data3file %>% filter(sampleNum == 2)

Barcodes = unique(S1$BC50StarcodeD8)

set.seed(2059)
subsample1 = sample(Barcodes,1000)
subsample2 = sample(Barcodes,1000)
subsample3 = sample(Barcodes,1000)
BarcodesLv1 = as.integer(stringdistmatrix(subsample1, method = "lv"))
BarcodesLv2 = as.integer(stringdistmatrix(subsample2, method = "lv"))
BarcodesLv3 = as.integer(stringdistmatrix(subsample3, method = "lv"))
lBarcodesLv = length(BarcodesLv1)

BarcodesLv = tibble(
  lvdist = c(BarcodesLv1, BarcodesLv2, BarcodesLv3),
  subsamNum = c(rep("subsamping1", lBarcodesLv), rep("subsamping2", lBarcodesLv), rep("subsamping3", lBarcodesLv)))

BarcodesLvHist <- BarcodesLv %>% group_by(subsamNum, lvdist) %>% summarise(length(lvdist)) %>%
  group_by(subsamNum) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

meanLvHisr = BarcodesLvHist %>% group_by(subsamNum) %>% summarise(mean = mean(lvdist)) #25, 24.5,24.5
ErrorMeanLvHisrProcessed50d8S1 = meanLvHisr %>% mutate(errorPercent = (100*(25-mean)/25), condition = "S150d8Processed")


FM01_ProcessedBarcodesLvHistPlot <- ggplot(BarcodesLvHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  facet_wrap(facets = vars(subsamNum)) +
  theme_classic() + 
  xlim(0,50) +
  geom_vline(aes(xintercept= mean),meanLvHisr, color="blue", size=0.5) +
  geom_vline(xintercept= 25, color="red", size=0.5)

FM01_ProcessedBarcodesLvHistPlotMagnified <- ggplot(BarcodesLvHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  facet_wrap(facets = vars(subsamNum)) +
  theme_classic() + 
  xlim(0,30) +
  ylim(0,0.03)

ggsave(FM01_ProcessedBarcodesLvHistPlot, file = paste0(plotDirectory, 'FM01_ProcessedBarcodesLvHistPlotS1.svg'), width = 10, height = 4)
ggsave(FM01_ProcessedBarcodesLvHistPlotMagnified, file = paste0(plotDirectory, 'FM01_ProcessedBarcodesLvHistPlotMagnifiedS1.svg'), width = 10, height = 4)

##S2 50d8
Barcodes = unique(S2$BC50StarcodeD8)

set.seed(2059)
subsample1 = sample(Barcodes,1000)
subsample2 = sample(Barcodes,1000)
subsample3 = sample(Barcodes,1000)
BarcodesLv1 = as.integer(stringdistmatrix(subsample1, method = "lv"))
BarcodesLv2 = as.integer(stringdistmatrix(subsample2, method = "lv"))
BarcodesLv3 = as.integer(stringdistmatrix(subsample3, method = "lv"))
lBarcodesLv = length(BarcodesLv1)

BarcodesLv = tibble(
  lvdist = c(BarcodesLv1, BarcodesLv2, BarcodesLv3),
  subsamNum = c(rep("subsamping1", lBarcodesLv), rep("subsamping2", lBarcodesLv), rep("subsamping3", lBarcodesLv)))

BarcodesLvHist <- BarcodesLv %>% group_by(subsamNum, lvdist) %>% summarise(length(lvdist)) %>%
  group_by(subsamNum) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

meanLvHisr = BarcodesLvHist %>% group_by(subsamNum) %>% summarise(mean = mean(lvdist)) #25, 24.5,24.5
ErrorMeanLvHisrProcessed50d8S2 = meanLvHisr %>% mutate(errorPercent = (100*(25-mean)/25), condition = "S250d8Processed")

###################################Same Analysis for BC40

data4file = as_tibble(read.table(paste0(macDirectory, input3Directory, 'stepFourStarcodeShavedReads40.txt'), stringsAsFactors=F, header = T))

S1 = data4file %>% filter(sampleNum == 1)
S2 = data4file %>% filter(sampleNum == 2)

Barcodes = unique(S1$BC40StarcodeD8)

set.seed(2059)
subsample1 = sample(Barcodes,1000)
subsample2 = sample(Barcodes,1000)
subsample3 = sample(Barcodes,1000)
BarcodesLv1 = as.integer(stringdistmatrix(subsample1, method = "lv"))
BarcodesLv2 = as.integer(stringdistmatrix(subsample2, method = "lv"))
BarcodesLv3 = as.integer(stringdistmatrix(subsample3, method = "lv"))
lBarcodesLv = length(BarcodesLv1)

BarcodesLv = tibble(
  lvdist = c(BarcodesLv1, BarcodesLv2, BarcodesLv3),
  subsamNum = c(rep("subsamping1", lBarcodesLv), rep("subsamping2", lBarcodesLv), rep("subsamping3", lBarcodesLv)))

BarcodesLvHist <- BarcodesLv %>% group_by(subsamNum, lvdist) %>% summarise(length(lvdist)) %>%
  group_by(subsamNum) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

meanLvHisr = BarcodesLvHist %>% group_by(subsamNum) %>% summarise(mean = mean(lvdist)) #21, 21.5,20.5
ErrorMeanLvHisrProcessed40d8S1 = meanLvHisr %>% mutate(errorPercent = (100*(20-mean)/20), condition = "S140d8Processed")


FM01_ProcessedBarcodesLvHistPlot <- ggplot(BarcodesLvHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  facet_wrap(facets = vars(subsamNum)) +
  theme_classic() + 
  xlim(0,40) +
  geom_vline(aes(xintercept= mean),meanLvHisr, color="blue", size=0.5) +
  geom_vline(xintercept= 20, color="red", size=0.5)

FM01_ProcessedBarcodesLvHistPlotMagnified <- ggplot(BarcodesLvHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  facet_wrap(facets = vars(subsamNum)) +
  theme_classic() + 
  xlim(0,25) +
  ylim(0,0.03)

ggsave(FM01_ProcessedBarcodesLvHistPlot, file = paste0(plotDirectory, 'FM01_ProcessedBarcodesLvHistPlot40.svg'), width = 10, height = 4)
ggsave(FM01_ProcessedBarcodesLvHistPlotMagnified, file = paste0(plotDirectory, 'FM01_ProcessedBarcodesLvHistPlotMagnified40.svg'), width = 10, height = 4)

Barcodes = unique(S2$BC40StarcodeD8)

set.seed(2059)
subsample1 = sample(Barcodes,1000)
subsample2 = sample(Barcodes,1000)
subsample3 = sample(Barcodes,1000)
BarcodesLv1 = as.integer(stringdistmatrix(subsample1, method = "lv"))
BarcodesLv2 = as.integer(stringdistmatrix(subsample2, method = "lv"))
BarcodesLv3 = as.integer(stringdistmatrix(subsample3, method = "lv"))
lBarcodesLv = length(BarcodesLv1)

BarcodesLv = tibble(
  lvdist = c(BarcodesLv1, BarcodesLv2, BarcodesLv3),
  subsamNum = c(rep("subsamping1", lBarcodesLv), rep("subsamping2", lBarcodesLv), rep("subsamping3", lBarcodesLv)))

BarcodesLvHist <- BarcodesLv %>% group_by(subsamNum, lvdist) %>% summarise(length(lvdist)) %>%
  group_by(subsamNum) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

meanLvHisr = BarcodesLvHist %>% group_by(subsamNum) %>% summarise(mean = mean(lvdist)) #21, 21.5,20.5
ErrorMeanLvHisrProcessed40d8S2 = meanLvHisr %>% mutate(errorPercent = (100*(20-mean)/20), condition = "S240d8Processed")
###################################Same Analysis for BC30

data5file = as_tibble(read.table(paste0(macDirectory, input3Directory, 'stepFourStarcodeShavedReads30.txt'), stringsAsFactors=F, header = T))

S1 = data5file %>% filter(sampleNum == 1)
S2 = data5file %>% filter(sampleNum == 2)

Barcodes = unique(S1$BC30StarcodeD8)

set.seed(2059)
subsample1 = sample(Barcodes,1000)
subsample2 = sample(Barcodes,1000)
subsample3 = sample(Barcodes,1000)
BarcodesLv1 = as.integer(stringdistmatrix(subsample1, method = "lv"))
BarcodesLv2 = as.integer(stringdistmatrix(subsample2, method = "lv"))
BarcodesLv3 = as.integer(stringdistmatrix(subsample3, method = "lv"))
lBarcodesLv = length(BarcodesLv1)

BarcodesLv = tibble(
  lvdist = c(BarcodesLv1, BarcodesLv2, BarcodesLv3),
  subsamNum = c(rep("subsamping1", lBarcodesLv), rep("subsamping2", lBarcodesLv), rep("subsamping3", lBarcodesLv)))

BarcodesLvHist <- BarcodesLv %>% group_by(subsamNum, lvdist) %>% summarise(length(lvdist)) %>%
  group_by(subsamNum) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

meanLvHisr = BarcodesLvHist %>% group_by(subsamNum) %>% summarise(mean = mean(lvdist)) #17.5, 17.5,17.5
ErrorMeanLvHisrProcessed30d8S1 = meanLvHisr %>% mutate(errorPercent = (100*(15-mean)/15), condition = "S130d8Processed")


FM01_ProcessedBarcodesLvHistPlot <- ggplot(BarcodesLvHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  facet_wrap(facets = vars(subsamNum)) +
  theme_classic() + 
  xlim(0,30) +
  geom_vline(aes(xintercept= mean),meanLvHisr, color="blue", size=0.5) +
  geom_vline(xintercept= 15, color="red", size=0.5)

FM01_ProcessedBarcodesLvHistPlotMagnified <- ggplot(BarcodesLvHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  facet_wrap(facets = vars(subsamNum)) +
  theme_classic() + 
  xlim(0,20) +
  ylim(0,0.03)

ggsave(FM01_ProcessedBarcodesLvHistPlot, file = paste0(plotDirectory, 'FM01_ProcessedBarcodesLvHistPlot30S1.svg'), width = 10, height = 4)
ggsave(FM01_ProcessedBarcodesLvHistPlotMagnified, file = paste0(plotDirectory, 'FM01_ProcessedBarcodesLvHistPlotMagnified30S1.svg'), width = 10, height = 4)

Barcodes = unique(S2$BC30StarcodeD8)

set.seed(2059)
subsample1 = sample(Barcodes,1000)
subsample2 = sample(Barcodes,1000)
subsample3 = sample(Barcodes,1000)
BarcodesLv1 = as.integer(stringdistmatrix(subsample1, method = "lv"))
BarcodesLv2 = as.integer(stringdistmatrix(subsample2, method = "lv"))
BarcodesLv3 = as.integer(stringdistmatrix(subsample3, method = "lv"))
lBarcodesLv = length(BarcodesLv1)

BarcodesLv = tibble(
  lvdist = c(BarcodesLv1, BarcodesLv2, BarcodesLv3),
  subsamNum = c(rep("subsamping1", lBarcodesLv), rep("subsamping2", lBarcodesLv), rep("subsamping3", lBarcodesLv)))

BarcodesLvHist <- BarcodesLv %>% group_by(subsamNum, lvdist) %>% summarise(length(lvdist)) %>%
  group_by(subsamNum) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

meanLvHisr = BarcodesLvHist %>% group_by(subsamNum) %>% summarise(mean = mean(lvdist)) #17.5, 17.5,17.5
ErrorMeanLvHisrProcessed30d8S2 = meanLvHisr %>% mutate(errorPercent = (100*(15-mean)/15), condition = "S230d8Processed")



# temp = data5file  %>% 
#   dplyr::rename(cellID = V1,
#                 UMI = V2,
#                 BC30StarcodeD8 = V6,
#                 sampleNum = V8) %>%
#   select(-c(V3,V4,V5,V7))
# 
# temp30 = temp %>%
#   unique() %>%
#   select(-UMI)
# 
# write.table(temp30, file= paste0(input3Directory,'stepFourStarcodeShavedReads30.txt'),row.names=F,col.names=T,quote=F,sep="\t")

#######

ErrorMeanLvHisr = bind_rows(ErrorMeanLvHisrS1,ErrorMeanLvHisrS2, ErrorMeanLvHisrProcessed50d8S1,ErrorMeanLvHisrProcessed50d8S2,ErrorMeanLvHisrProcessed40d8S1, ErrorMeanLvHisrProcessed40d8S2, ErrorMeanLvHisrProcessed30d8S1, ErrorMeanLvHisrProcessed30d8S2)
ErrorMeanLvHisr$condition <- factor(ErrorMeanLvHisr$condition, levels=c("S1Raw", "S2Raw", "S150d8Processed", "S250d8Processed", "S140d8Processed", "S240d8Processed", "S130d8Processed", "S230d8Processed"))

tableAllMeanSEM = ErrorMeanLvHisr %>% group_by(condition) %>% summarise(mean = mean(errorPercent), sem = sqrt(var(errorPercent)/length(errorPercent))) %>% mutate(name = 1)

write.table(ErrorMeanLvHisr, file=paste0(plotDirectory,'FM01_ErrorMeanLvHisr.tsv'), col.names = TRUE, sep='\t')


plotErrorAll <- ggplot() +
  geom_jitter(aes(x= ErrorMeanLvHisr$condition,y= ErrorMeanLvHisr$errorPercent),alpha=0.6, size = 2.5, shape = 16, width = 0.1) +
  geom_pointrange(aes(x= tableAllMeanSEM$condition, y = tableAllMeanSEM$mean, ymin=tableAllMeanSEM$mean-tableAllMeanSEM$sem, ymax=tableAllMeanSEM$mean +tableAllMeanSEM$sem),size = 0.8, shape = 16)+
  ylim(-25,25) +
  geom_hline(aes(yintercept= 0)) +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x=element_blank())

ggsave(plotErrorAll, file = paste0(plotDirectory, 'FM01_plotErrorAll.svg'), width = 8, height = 4)

 
