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


homeDirectory <- '~/Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM01/'
plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/analysisPlots/"

barcode50 = as_tibble(read.table(paste0(homeDirectory, 'stepFourStarcodeShavedReads50.txt'), stringsAsFactors=F, header = T))
umiCut = 4 # Minimum UMI cutoff for reliable analysis. Can be played around with between 3 to 5.
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.  

upToLineageCounting = barcode50 %>% 
  filter(sampleNum %in% c("1", "2")) %>%
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
  filter(nLineages <= linCut) %>% select(-nUMI, -nLineages) %>%
  unique()


linCountTooverlapsColony = linCountTooverlaps %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))

ColonyDistS1 = linCountTooverlapsColony %>% filter(sampleNum == "S1")
ColonyDistS2 = linCountTooverlapsColony %>% filter(sampleNum == "S2")

cloneDistribution = ggplot(ColonyDistS1, aes(nColony)) +
  geom_histogram(aes(x=nColony), binwidth = 1,color="black", fill="black",boundary = 0)+
  geom_rug() + 
  theme_classic()+
  theme(text = element_text(size = 22))+ 
  ylab("Frequency")+
  xlab('Clone Size') + 
  scale_x_continuous(breaks = seq(0,160,40))+
  annotate("text", x=80, y=2000, label= paste0("Total barcoded cells = ", sum(ColonyDistS1$nColony)), size = 6) + 
  annotate("text", x=80, y=1930, label= paste0("Total number of singletons = ", sum(ColonyDistS1$nColony==1), "  (",round(100*length(ColonyDistS1$nColony==1)/ sum(ColonyDistS1$nColony),2),"%)"), size=6)+
  annotate("text", x=80, y=1860, label= paste0("Largest clone fraction = ", round(max(ColonyDistS1$nColony)/sum(ColonyDistS1$nColony),4) ), size=6)
ggsave(cloneDistribution, file = paste0(plotDirectory, 'FM01_CloneSizeHistogramS1.svg'), width = 10, height = 10)

cloneDistribution2 = ggplot(ColonyDistS2, aes(nColony)) +
  geom_histogram(aes(x=nColony), binwidth = 1,color="black", fill="black",boundary = 0)+
  geom_rug() + 
  theme_classic()+
  theme(text = element_text(size = 22))+ 
  ylab("Frequency")+
  xlab('Clone Size') + 
  scale_x_continuous(breaks = seq(0,160,40))+
  annotate("text", x=80, y=2000, label= paste0("Total barcoded cells = ", sum(ColonyDistS2$nColony)), size = 6) + 
  annotate("text", x=80, y=1930, label= paste0("Total number of singletons = ", sum(ColonyDistS2$nColony==1), "  (",round(100*length(ColonyDistS2$nColony==1)/ sum(ColonyDistS2$nColony),2),"%)"), size=6)+
  annotate("text", x=80, y=1860, label= paste0("Largest clone fraction = ", round(max(ColonyDistS2$nColony)/sum(ColonyDistS2$nColony),4) ), size=6)
ggsave(cloneDistribution2, file = paste0(plotDirectory, 'FM01_CloneSizeHistogramS2.svg'), width = 10, height = 10)
