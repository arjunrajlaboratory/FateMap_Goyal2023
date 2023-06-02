# make spheroid plots
library(dplyr)
library(ggplot2)
home1Directory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/microscopy/spheroids/'
plotDirectory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/plotData/analysisPlots/'

data2file = paste0(home1Directory, '20190906_Spheroid_QuantPlotAll.csv')

data2 = as_tibble(read.csv(data2file, stringsAsFactors=F, header = T))
dataSelect = data2 %>% filter(sampleName %in% c("Rc7_5nMT", "Rc11_5nMT", "WM989", "Rc21_5nMT"))
dataTram = data2 %>% filter(condition %in% c("5nMT", "naive"))
dataPLX = data2 %>% filter(condition %in% c("250nmPLX","1uMPLX", "naive"))

Spheroid <- ggplot(data = dataSelect, aes(x=reorder(sampleName, ratio, FUN = mean),y=ratio)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 2, shape = 16) +
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica"))+
  ylim(0, 1.1*(max(dataSelect$ratio)))
ggsave(Spheroid, file = paste0(plotDirectory, "SpheroidQuantification.svg"), width = 3, height = 4)

Spheroid <- ggplot(data = dataTram, aes(x=reorder(sampleName, ratio, FUN = mean),y=ratio)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1, shape = 16) +
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica"))+
  ylim(0, 1.05*(max(dataTram$ratio)))
ggsave(Spheroid, file = paste0(plotDirectory, "SpheroidQuantificationTramAll.svg"), width = 5, height = 4)

Spheroid <- ggplot(data = dataPLX, aes(x=reorder(sampleName, ratio, FUN = mean),y=ratio)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1, shape = 16) +
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica"))+
  ylim(0, 1.05*(max(dataPLX$ratio)))
ggsave(Spheroid, file = paste0(plotDirectory, "SpheroidQuantificationPLXAll.svg"), width = 7, height = 4)