########################################################################################
#Description: Analysis of clusters from DEGs for data from Herlyn Lab and Tirosh dataset
#Author: Maalavika Pillai
#Version: 1
#Edited on:14th December
########################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(svglite) 
library(ggrepel)
library(SeuratObject)
library(hdf5r)
library(VennDiagram)
library(RColorBrewer)
library(tidyr)

dataDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/MHdata/seuratOutput/"
dataDirectory1 <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/TiroshData/seuratOutput/"
plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/plots/"

MH_Data <- read.csv(paste0(dataDirectory, "MH_Random_real_overlap.csv"))
IT_Data <- read.csv(paste0(dataDirectory1, "IT_Random_real_overlap.csv"))
compiled_data <- rbind(data.frame(MH_Data, "Sample" = "MH"), data.frame(IT_Data, "Sample" = "IT"))
compiled_data_DEG <- compiled_data %>% 
  filter(name!="VarGenes")
mean_line <-compiled_data_DEG %>% group_by(Sample,sample) %>% summarise(mean_y = mean(value));
plot <- ggplot(compiled_data_DEG, aes(x = sample, y = value, col = sample))+
  geom_jitter( position=position_jitterdodge(dodge.width =0.8, jitter.width = 0.2), shape = 16) +
  labs(x="", y = "Fraction overlap")+
  facet_wrap(.~Sample)+
  geom_path(data = mean_line, aes(x = sample, y = mean_y, group = 1), color = "grey")+
  theme_classic() 
ggsave(plot = plot,filename= paste0(plotDirectory, "MH_IT_Random_real_overlap.svg"), height=10, width=10)
