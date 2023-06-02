#####################NOTE#####################
#####Plotting data collected by YG + GB
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/YGGB Paper Figures/FinalTables/scriptToPlot.R
##############################################

library(ggplot2)
library(tidyverse)
library(dplyr)

dataDirectory <- "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/microscopy/DAPIScans/WM989/drugDoseType/"
data2Directory <- "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/microscopy/DAPIScans/WM989/DOT1L/"
plotDirectory <- "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/plotData/analysisPlots/"
####First Low vs High Dose and High vs Trametenib
LowHighTram = as_tibble(read.table(file = paste0(dataDirectory, "0902_lowHighTram.csv"), header = TRUE, stringsAsFactors=F, sep = ","))

foldChangeLowHigh = LowHighTram %>% filter(condition %in% c("lowPLX", "highPLX")) %>% dplyr::select(Var1, totalcolonies_well, condition) %>% 
  spread(condition, totalcolonies_well) %>% 
  mutate(colonyFoldChange = lowPLX/highPLX) %>% mutate(name = 1)

tableAllMeanSEM = foldChangeLowHigh %>% summarise(mean = mean(colonyFoldChange), sem = sqrt(var(colonyFoldChange)/length(colonyFoldChange))) %>% mutate(name = 1)

plot <- ggplot() +
  geom_jitter(aes(x= foldChangeLowHigh$name,y= foldChangeLowHigh$colonyFoldChange),alpha=0.6, size = 2.5, shape = 16, width = 0.1) +
  geom_pointrange(aes(x= tableAllMeanSEM$name, y = tableAllMeanSEM$mean, ymin=tableAllMeanSEM$mean-tableAllMeanSEM$sem, ymax=tableAllMeanSEM$mean +tableAllMeanSEM$sem),size = 0.8, shape = 16)+
  ylim(0,1.1*max(foldChangeLowHigh$colonyFoldChange))+
  xlim(0.8,1.2) +
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(plot = plot, file = paste0(plotDirectory, 'foldChangeColonyDose','.svg'), width = 1.5, height = 4)

####Trametenib Plot
foldChangeHighTramColony = LowHighTram %>% filter(condition %in% c("Tram", "highPLX")) %>% dplyr::select(Var1, totalcolonies_well, condition) %>% 
  spread(condition, totalcolonies_well) %>% 
  mutate(colonyFoldChange = highPLX/Tram) %>% mutate(name = 1)

foldChangeHighTramCells = LowHighTram %>% filter(condition %in% c("Tram", "highPLX")) %>% dplyr::select(Var1, cells_outside, condition) %>% 
  spread(condition, cells_outside) %>% 
  mutate(colonyFoldChange = highPLX/Tram) %>% mutate(name = 2)

foldChangeHighTram = bind_rows(foldChangeHighTramColony, foldChangeHighTramCells)

tableAllMeanSEM = foldChangeHighTram %>% group_by(name) %>% summarise(mean = mean(colonyFoldChange), sem = sqrt(var(colonyFoldChange)/length(colonyFoldChange))) 

plot <- ggplot() +
  geom_jitter(aes(x= foldChangeHighTram$name,y= foldChangeHighTram$colonyFoldChange),alpha=0.6, size = 2.5, shape = 16, width = 0.1) +
  geom_pointrange(aes(x= tableAllMeanSEM$name, y = tableAllMeanSEM$mean, ymin=tableAllMeanSEM$mean-tableAllMeanSEM$sem, ymax=tableAllMeanSEM$mean +tableAllMeanSEM$sem),size = 0.8, shape = 16)+
  ylim(0,1.1*max(foldChangeHighTram$colonyFoldChange))+
  xlim(0.8,2.2) +
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(plot = plot, file = paste0(plotDirectory, 'foldChangeColonyCellsDrugType','.svg'), width = 3, height = 4)

####DOT1L vs DMSO
DOT1LCtrl = as_tibble(read.table(file = paste0(data2Directory, "WellA3_B1FinalDataTable.csv"), header = TRUE, stringsAsFactors=F, sep = ","))

foldChangeCtrlDot1LColony = DOT1LCtrl %>% filter(condition %in% c("Ctrl", "DOT1L")) %>% dplyr::select(Var1, totalcolonies_well, condition) %>% 
  spread(condition, totalcolonies_well) %>% 
  mutate(colonyFoldChange = DOT1L/Ctrl) %>% mutate(name = 1)

foldChangeCtrlDot1LCells = DOT1LCtrl %>% filter(condition %in% c("Ctrl", "DOT1L")) %>% dplyr::select(Var1, cells_outside, condition) %>% 
  spread(condition, cells_outside) %>% 
  mutate(colonyFoldChange = DOT1L/Ctrl) %>% mutate(name = 2)

foldChangeCtrlDot1 = bind_rows(foldChangeCtrlDot1LColony, foldChangeCtrlDot1LCells)

tableAllMeanSEM = foldChangeCtrlDot1 %>% group_by(name) %>% summarise(mean = mean(colonyFoldChange), sem = sqrt(var(colonyFoldChange)/length(colonyFoldChange))) 

plot <- ggplot() +
  geom_jitter(aes(x= foldChangeCtrlDot1$name,y= foldChangeCtrlDot1$colonyFoldChange),alpha=0.6, size = 2.5, shape = 16, width = 0.1) +
  geom_pointrange(aes(x= tableAllMeanSEM$name, y = tableAllMeanSEM$mean, ymin=tableAllMeanSEM$mean-tableAllMeanSEM$sem, ymax=tableAllMeanSEM$mean +tableAllMeanSEM$sem),size = 0.8, shape = 16)+
  ylim(0,1.1*max(foldChangeCtrlDot1$colonyFoldChange))+
  xlim(0.8,2.2) +
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(plot = plot, file = paste0(plotDirectory, 'foldChangeColonyCellsDOT1LCtrl','.svg'), width = 3, height = 4)
