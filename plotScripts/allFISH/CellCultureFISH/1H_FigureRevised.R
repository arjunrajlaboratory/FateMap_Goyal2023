########################################################################################
#Description: Spot counts Representation
#Author: YG
#Version: 1
#Created on: February 2 2023
########################################################################################


library(dplyr)
library(tidyr)
library(svglite)
library(FactoMineR)
library(ggplot2)
library(ggpubr)

plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/1H/Plots/New_Tito/"

plot1Directory <- "/Users/yogeshgoyal/Dropbox (RajLab)/FateMap/REVISIONS/analysisPlots/"

avg_geneCounts = as_tibble(read.table(file = paste0(plotDirectory, "AverageGeneCounts.csv"), header = TRUE, stringsAsFactors=F, sep = ","))

avg_geneCounts1 = avg_geneCounts %>% select(-YFP) %>% gather(key = "gene" , value = "count", 2:4)





###
plot = ggplot() +
  geom_point(data = avg_geneCounts, aes(x = "gene", y = BGN), color = "gray50", size = 1.5, shape = 16) +
  theme_classic()
ggsave(plot, file = paste0(plot1Directory, '1H_BGN.svg'), width = 1, height = 3)

  
  
  
geneCounts.pca <- PCA(avg_geneCounts[2:4])
p <- plot.PCA(geneCounts.pca,label = "none")
p <-  p + theme_classic() + 
  border()
print(p)
ggsave(plot = p, filename = paste0(plotDirectory,"AvgSpotCount_PCA_all.svg"))
write.csv(avg_geneCounts, paste0(plotDirectory,"AverageGeneCounts.csv"))