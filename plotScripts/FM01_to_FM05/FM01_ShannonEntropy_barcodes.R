########################################################################################
#Description: Shannon entropy to check transcriptional homogeneity within clones
#Author: Maalavika Pillai
#Version: 1
#Edited on: 7th December 2022
########################################################################################
library(tidyverse)
library(ggplot2)
library(vegan)
library(ggpubr)
library(svglite)

dataDirectory <-"~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM01/Data/"
plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/Entropy/"

#Load data and barcodes
logNormalizedCounts = as_tibble(read.table(file = paste0(dataDirectory, "logNormalizedSCTCountsS1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
logNormalizedCounts_S1 = logNormalizedCounts %>% 
  filter(sampleNum=="S1") %>%
  select(-sampleNum)
barcode50 = as_tibble(read.table(paste0(dataDirectory, 'barcodeCellID.tsv'), stringsAsFactors=F, header = T))
barcodes_S1 = barcode50 %>% 
  filter(sampleNum=="S1"&cellID %in% logNormalizedCounts_S1$cellID) %>%
  select(-sampleNum) %>%
  group_by(BC50StarcodeD8) %>%
  mutate("nSize" = length(cellID)) %>%
  filter(nSize > 4) %>%
  select(-nSize)
cellIDs <- barcodes_S1$cellID #Only selects barcoded cells

clusters_sub <- read.delim(paste0(dataDirectory, "umapClusters_s1s2Scanorama_50pcs_filter_snn0_4.tsv"), sep = "\t")
clusters_sub <- clusters_sub %>% filter(sampleNum=="S1")
rownames(clusters_sub) <- clusters_sub$cellID
clusters <- clusters_sub[-c(3)]
for (res in c(6,8)) {
  clusters_sub <- read.delim(paste0(dataDirectory, "umapClusters_s1s2Scanorama_50pcs_filter_snn0_",res,".tsv"), sep = "\t")
  clusters_sub <- clusters_sub %>% filter(sampleNum=="S1")
  rownames(clusters_sub) <- clusters_sub$cellID
  clusters_sub <- clusters_sub[rownames(clusters),-c(2:3)]
  clusters <- cbind(clusters, clusters_sub)
}
names(clusters)[c(1,3,4)] <- c("seurat_cluster_0.4","seurat_cluster_0.6","seurat_cluster_0.8")
clusters <- clusters[cellIDs,]
logNormalizedCounts_S1 <- logNormalizedCounts_S1 %>% filter(cellID %in% cellIDs)

entropy_combined <- data.frame()
for(cluster_res in c(1,3,4)){
  seurat_clusters <- clusters[,c("cellID", names(clusters)[cluster_res])]
  names(seurat_clusters) <- c("cellID", "seurat_clusters")
  clus_geneexp_barcode <- inner_join(logNormalizedCounts_S1, seurat_clusters, by = "cellID")
  clus_geneexp_barcode <- inner_join(clus_geneexp_barcode, barcodes_S1)
  
  
  #Entropy for each cluster
  entropy_clones <- clus_geneexp_barcode %>%
    group_by(BC50StarcodeD8) %>%
    count(seurat_clusters) %>%
    pivot_wider(names_from = seurat_clusters, values_from = n)
  entropy_clones[is.na(entropy_clones)] <- 0
  entropy_clones_val <- apply(entropy_clones, 1, function(x){diversity(as.numeric(x[-1]))})
  
  
  
  #Shuffle barcodes for random assignment 
  set.seed(100)
  rand_clus_geneexp_barcode <- clus_geneexp_barcode %>%
    mutate(BC50StarcodeD8 = sample(BC50StarcodeD8))
  entropy_random <- rand_clus_geneexp_barcode %>%
    group_by(BC50StarcodeD8) %>%
    count(seurat_clusters) %>%
    pivot_wider(names_from = seurat_clusters, values_from = n)
  entropy_random[is.na(entropy_random)] <- 0
  entropy_random_val <- apply(entropy_random, 1, function(x){diversity(as.numeric(x[-1]))})
  # entropy_clones_val <- entropy_clones_val/log(ncol(entropy_clones)-1) #Shannon Equitability Index/ Normalized to 1
  # entropy_random_val <- entropy_random_val/ log(ncol(entropy_clones)-1)
  entropy_all <- data.frame("Clones" = entropy_clones_val,"Random" = entropy_random_val, "Barcodes" = entropy_clones$BC50StarcodeD8, 
                            "resolution" = gsub("seurat_cluster_","", names(clusters)[cluster_res]))
  entropy_all  <- entropy_all %>% pivot_longer(names_to = "Sample", cols = Clones:Random)
  entropy_combined <- rbind(entropy_combined, entropy_all)
}
entropy_combined$resolution <- as.factor(as.character(entropy_combined$resolution))
#Dataframe
#confirm order of barcodes is same/paired:
#all(entropy_random$BC50StarcodeD8==entropy_clones$BC50StarcodeD8)
plot <- ggplot(entropy_combined, aes(y = value, x = Sample, col = Sample)) + 
  geom_boxplot(width = 0.3) +
  geom_jitter(width = 0.1, shape = 16) +
  stat_compare_means( size =4, label.x.npc = "left", paired = T)+
  facet_wrap(.~resolution)+
  labs(x="", y = "Shannon Diversity Index")+
  theme_classic() + 
  theme(legend.position="none") 
ggsave(plot = plot,filename=paste0(plotDirectory, "ShannonDiversityIndex_resolutions.svg"), width= 18, height =12)

##############
#For normalized entropy - Shannon's Equitability Index

plot <-ggplot(entropy_combined, aes(y = value/log(ncol(entropy_clones)-1), x = Sample, col = Sample)) + 
  geom_boxplot(width = 0.3) +
  geom_jitter(width = 0.1, shape = 16) +
  stat_compare_means( size =4, label.x.npc = "left", paired = T)+
  facet_wrap(.~resolution)+
  labs(x="", y = "Shannon Diversity Index")+
  theme_classic() + 
  theme(legend.position="none") 
ggsave(plot = plot,filename=paste0(plotDirectory, "ShannonEquitabilityIndex.svg"), width= 12, height = 12)

