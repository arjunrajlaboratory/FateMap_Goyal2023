#####################NOTE#####################
#Code to check for E-distances and Knn clusters 
#First started by Jonas Braun on April 9, 2023; Last updated by Yogesh Goyal on April 10, 2023
##############################################

#Lineplot stats
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

dataDirectory <- "~/Dropbox (RajLab)/FateMap/FINAL_REVISION/analysis_Jonas/data/e_distResults/"
plotDirectory <- "~/Dropbox (RajLab)/FateMap/FINAL_REVISION/analysis_Jonas/plots/"

dfNumberClusters <- tibble(read.table(file = paste0(dataDirectory, "NumberClusters_fm01VSfm06.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))
dfDistance<- tibble(read.table(file = paste0(dataDirectory,"Edistances_fm01VSfm06.csv"), header = TRUE, stringsAsFactors=F, sep = ",", fill = TRUE))

dfMeanCluster <- dfNumberClusters %>% group_by(Dataset,resolution) %>% dplyr::summarise(mean = mean(numberCluster)) 
dfSDCluster <- dfNumberClusters %>% group_by(Dataset,resolution) %>% dplyr::summarise(sd = sd(numberCluster)) 
df_joinCluster <- dfMeanCluster %>% 
  left_join(dfSDCluster)

df_joinClusterSpread = df_joinCluster %>% dplyr::select(-sd) %>% spread(Dataset, mean)

wilcox.test(df_joinClusterSpread$fm01, df_joinClusterSpread$fm06, paired = TRUE, alternative = "two.sided", )


rm(dfMeanCluster)
rm(dfSDCluster)


dfMeanDistance <- dfDistance %>% group_by(Dataset,numberCluster) %>% dplyr::summarise(mean = mean(distances)) 
dfSDDistance <- dfDistance %>% group_by(Dataset, numberCluster) %>% dplyr::summarise(sd = sd(distances)) 
df_joinDistance <- dfMeanDistance %>% 
  left_join(dfSDDistance)



rm(dfMeanDistance)
rm(dfSDDistance)

#resolution <-> number of clusters
plot <- ggplot(data = dfNumberClusters, aes(x=resolution, y=numberCluster, color = Dataset)) +
  geom_jitter(height = 0.1, shape = 16) + 
  geom_line(data = df_joinCluster, mapping = aes(x=resolution, y=mean)) + 
  scale_x_continuous(limits = c(0,1.05), breaks = c(0.2, 0.4, 0.6, 0.8, 1)) +
  ylab("number of clusters") + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave(plot, file = paste0(plotDirectory, "resolution_NumberCluster_FM01FM06.svg"), width = 6, height = 4)

# #number of cluster <-> e-dist; all points
# plot <- ggplot(data = dfDistance, aes(x=numberCluster, y=distances, color = Dataset)) +
#   geom_point(position = position_dodge(0.5)) + 
#   geom_line(data = df_joinDistance, mapping = aes(x=numberCluster, y=mean)) + 
#   scale_x_continuous(limits = c(4,16)) +
#   ylab("E-distance") + xlab("number of clusters") + 
#   theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
#   theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
#         panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
# ggsave(plot, file = paste0(plotDirectory, "numberCluster_Edist_all_FM01FM06.svg"), width = 10, height = 6)




#connecting dots
dfDistanceSample<- dfDistance %>% filter(sample ==1)

dfDistancePlot <- data.frame()
for (i in c(6,8,10,12,14))
{
  dfDistanceFilter1 <- dfDistanceSample %>% filter((Dataset== "fm01") & (numberCluster == i))
  dfDistanceFilter6 <- dfDistanceSample %>% filter((Dataset== "fm06") & (numberCluster == i))
  val1 <- as.array(unique(dfDistanceFilter1$id))
  val6 <- as.array(unique(dfDistanceFilter6$id))
  nSimul <- min(length(unique(val1)), length(unique(val6)))
  
  if (length(val1) <= 1)
  {
    id1 = val1
  }else
  {
    id1 <- sample(val1, nSimul)
  }
  
  if (length(val6) <= 1)
  {
    id6 <- val6
  }else
  {
    id6 <- sample(val6, nSimul)
  }
  dfDistanceFilter1 <- dfDistanceFilter1 %>% filter(id %in% id1)
  dfDistanceFilter6 <- dfDistanceFilter6 %>% filter(id %in% id6)
  dfDistanceFilter1$rank <- rank(-dfDistanceFilter1$distances, ties.method = "min")
  dfDistanceFilter6$rank <- rank(-dfDistanceFilter6$distances, ties.method = "min")
  dfDistancePlot <- do.call("rbind", list(dfDistancePlot, dfDistanceFilter1, dfDistanceFilter6))
}
rm(dfDistanceFilter1)
rm(dfDistanceFilter6)
rm(dfDistanceSample)
rm(id1)
rm(id6)
rm(nSimul)
rm(i)
rm(val1)
rm(val6)

dfDistancePlot$clusterRank = paste(dfDistancePlot$numberCluster, dfDistancePlot$rank, sep="_")


plot = ggplot(dfDistancePlot,aes(x=Dataset, y=distances, color = Dataset)) +
  geom_point(position = position_dodge(0.5)) +
  geom_line(aes(group=rank), alpha = 0.2, color = "black") +
  facet_wrap(~numberCluster, nrow = 1) +
  stat_compare_means(paired = TRUE,method.args = list(alternative = "less"), tip.length = 0.02, size = 3)+
  scale_color_manual(values = c(blues9[8], blues9[4]))+
  theme_classic() +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plot, file = paste0(plotDirectory, "numberCluster_Edist_connectingLines_FM01FM06.svg"), width = 6, height = 4)