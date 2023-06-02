#####################NOTE#####################
#Cleaned up version of V1:
##############################################

######Things to do tomorrow
#Create the manually curated set of randomized-morphologyMap (/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/bulkRNASeq/RandomizationRevised.tsv)
##create a new file that merges this ^^ file and (/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/bulkRNASeq/BulkColonySeq_20200220_meltedDataAll_normalized_hg19.tsv)
##run the following script. remove ones that are uncategorized. 

####script to analyze/heatmap for fast vs slow invadig cells based on individual colonies
library(dplyr)
library(ggplot2)
library(gplots)
library(tidyverse)
library(RColorBrewer)

dataDirectory <- "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/bulkRNASeq/"
home1Directory <-("/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM01/")
plotDirectory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/plotData/analysisPlots/'

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

######combining Data
RNASeqCounts = as_tibble(read.table(file = paste0(dataDirectory,"BulkColonySeq_20200220_meltedDataAll_normalized_hg19.tsv"), header = TRUE, sep = "\t"))
Random_demultiplex_Revised = as_tibble(read.table(file = paste0(dataDirectory,"RandomizationRevised.tsv"), header = TRUE, sep = "\t"))
Random_demultiplex_Revised$Original = sub("^", "sample_", Random_demultiplex_Revised$Original)

RNASeqCounts1 = RNASeqCounts %>% mutate(sampleID = gsub("WM989-ResistantColonies-","",RNASeqCounts$sampleID),
                                        sampleID = as.integer(sampleID)) 
RNASeqCountsRevised = inner_join(RNASeqCounts1,Random_demultiplex_Revised, by = "sampleID") %>% dplyr::select(-sampleID) %>% arrange(Original)
write_tsv(RNASeqCountsRevised, paste0(dataDirectory,"deRandomized_BulkColonySeq_20200220_meltedData_hg19_normalizedRevised.tsv"), col_names = T)

#counts table for all XX samples
counts = as_tibble(read.table(file = paste0(dataDirectory,"deRandomized_BulkColonySeq_20200220_meltedData_hg19_normalizedRevised.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
counts1 = counts %>% dplyr::select(tpm, gene_name,Original, colonyTypeFine) %>% filter(!colonyTypeFine == "naive") %>% filter(!colonyTypeFine == "uncategorized") 
counts1$colonyTypeFine = recode(counts1$colonyTypeFine, notOnTop = "Rest")
counts1$colonyTypeFine = recode(counts1$colonyTypeFine, onTop = "Rest")

countsFilter1 <- aggregate(counts1['tpm'], by=counts1[c('Original', 'gene_name', 'colonyTypeFine')], sum)

avgTPM = countsFilter1 %>% group_by(colonyTypeFine, gene_name) %>% summarise(avgTPM = mean(tpm))

allCounts_rawCountsMatrixNew = avgTPM %>% spread(colonyTypeFine, avgTPM)
allCounts_rawCountsMatrixNew[is.na(allCounts_rawCountsMatrixNew)] <- 0

foldChangeAvgTPM = allCounts_rawCountsMatrixNew %>% filter(Rest+small >10) %>% 
  mutate(log2FC = log2((small+1)/(Rest+1))) %>%
  filter(log2FC >1.5| log2FC < -1.5)

s1s2_scTransform.markers = as_tibble(read.table(file = paste0(plotDirectory, "scTransformMarkers_snn06.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
s1s2_scTransform.markersSubset <- s1s2_scTransform.markers %>% filter(p_val_adj <0.05, avg_logFC>1)

geneList10X = as_tibble(unique(s1s2_scTransform.markersSubset$gene))
geneListBulk = as_tibble(unique(foldChangeAvgTPM$gene_name))
s1s2_scTransform.markersSubsetBulkSubset = s1s2_scTransform.markersSubset %>% filter(gene %in% geneListBulk$value)

avgTPMselect = foldChangeAvgTPM %>% filter(gene_name %in% geneList10X$value)
upRegulatedGeneList = avgTPMselect %>% filter(log2FC>1.5) %>% select(gene_name)
downRegulatedGeneList = avgTPMselect %>% filter(log2FC< -1.5) %>% select(gene_name)

testCoverageA = s1s2_scTransform.markersSubset %>% group_by(cluster) %>% summarise(fraction = (sum(gene %in% upRegulatedGeneList$gene_name))/(length(gene))) %>% mutate(normalizedFraction = fraction/(max(fraction)))
testCoverage2A = s1s2_scTransform.markersSubset %>% group_by(cluster) %>% summarise(fraction = (sum(gene %in% downRegulatedGeneList$gene_name))/(length(gene))) %>% mutate(normalizedFraction = fraction/(max(fraction)))

umapCoordinates = as_tibble(read.table(file = paste0(home1Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapClusters = as_tibble(read.table(file = paste0(home1Directory, "umapClusters_s1s2Scanorama_50pcs_filter_snn0_6.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% rename(cluster = seurat_clusters)

umapClustersFractionUp = inner_join(umapClusters,testCoverageA, by = "cluster")

plot <- ggplot(umapCoordinates, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes_string(color = umapClustersFractionUp$normalizedFraction), size = 1, shape = 16) +
  scale_color_gradient(low = "gray93", high = "red4") +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Scaled\n Log UMI counts")
ggsave(plot, file = paste0(plotDirectory, 'phenotypeNGFRColonyMap.svg'), width = 6, height = 6.429)

write.table(avgTPMselect, file=paste0(plotDirectory,'avgTPMselectSmall.tsv'), col.names = TRUE, sep='\t')

####################################
##ON Top
counts1 = counts %>% dplyr::select(tpm, gene_name,Original, colonyTypeFine) %>% filter(!colonyTypeFine == "naive") %>% filter(!colonyTypeFine == "uncategorized") 
counts1$colonyTypeFine = recode(counts1$colonyTypeFine, notOnTop = "Rest")
counts1$colonyTypeFine = recode(counts1$colonyTypeFine, small = "Rest")

countsFilter1 <- aggregate(counts1['tpm'], by=counts1[c('Original', 'gene_name', 'colonyTypeFine')], sum)

avgTPM = countsFilter1 %>% group_by(colonyTypeFine, gene_name) %>% summarise(avgTPM = mean(tpm))

allCounts_rawCountsMatrixNew = avgTPM %>% spread(colonyTypeFine, avgTPM)
allCounts_rawCountsMatrixNew[is.na(allCounts_rawCountsMatrixNew)] <- 0

foldChangeAvgTPM = allCounts_rawCountsMatrixNew %>% filter(Rest+onTop >10) %>% 
  mutate(log2FC = log2((onTop+1)/(Rest+1))) %>%
  filter(log2FC >1.5| log2FC < -1.5)

s1s2_scTransform.markers = as_tibble(read.table(file = paste0(plotDirectory, "scTransformMarkers_snn06.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
s1s2_scTransform.markersSubset <- s1s2_scTransform.markers %>% filter(p_val_adj <0.05, avg_logFC>1)

geneList10X = as_tibble(unique(s1s2_scTransform.markersSubset$gene))
geneListBulk = as_tibble(unique(foldChangeAvgTPM$gene_name))
s1s2_scTransform.markersSubsetBulkSubset = s1s2_scTransform.markersSubset %>% filter(gene %in% geneListBulk$value)

avgTPMselect = foldChangeAvgTPM %>% filter(gene_name %in% geneList10X$value)
upRegulatedGeneList = avgTPMselect %>% filter(log2FC>1.5) %>% select(gene_name)
downRegulatedGeneList = avgTPMselect %>% filter(log2FC< -1.5) %>% select(gene_name)

testCoverageA = s1s2_scTransform.markersSubset %>% group_by(cluster) %>% summarise(fraction = (sum(gene %in% upRegulatedGeneList$gene_name))/(length(gene))) %>% mutate(normalizedFraction = fraction/(max(fraction)))
testCoverage2A = s1s2_scTransform.markersSubset %>% group_by(cluster) %>% summarise(fraction = (sum(gene %in% downRegulatedGeneList$gene_name))/(length(gene))) %>% mutate(normalizedFraction = fraction/(max(fraction)))

umapCoordinates = as_tibble(read.table(file = paste0(home1Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapClusters = as_tibble(read.table(file = paste0(home1Directory, "umapClusters_s1s2Scanorama_50pcs_filter_snn0_6.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% rename(cluster = seurat_clusters)

umapClustersFractionUp = inner_join(umapClusters,testCoverageA, by = "cluster")

plot <- ggplot(umapCoordinates, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes_string(color = umapClustersFractionUp$normalizedFraction), size = 1, shape = 16) +
  scale_color_gradient(low = "gray93", high = "red4") +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Scaled\n Log UMI counts")
ggsave(plot, file = paste0(plotDirectory, 'phenotypeVCAM1ColonyMap.svg'), width = 6, height = 6.429)

write.table(avgTPMselect, file=paste0(plotDirectory,'avgTPMselectOnTop.tsv'), col.names = TRUE, sep='\t')

######
counts1 = counts %>% dplyr::select(tpm, gene_name,Original, colonyTypeFine) %>% filter(!colonyTypeFine == "naive") %>% filter(!colonyTypeFine == "uncategorized") 
counts1$colonyTypeFine = recode(counts1$colonyTypeFine, onTop = "Rest")
counts1$colonyTypeFine = recode(counts1$colonyTypeFine, small = "Rest")

countsFilter1 <- aggregate(counts1['tpm'], by=counts1[c('Original', 'gene_name', 'colonyTypeFine')], sum)

avgTPM = countsFilter1 %>% group_by(colonyTypeFine, gene_name) %>% summarise(avgTPM = mean(tpm))

allCounts_rawCountsMatrixNew = avgTPM %>% spread(colonyTypeFine, avgTPM)
allCounts_rawCountsMatrixNew[is.na(allCounts_rawCountsMatrixNew)] <- 0

foldChangeAvgTPM = allCounts_rawCountsMatrixNew %>% filter(Rest+notOnTop >10) %>% 
  mutate(log2FC = log2((notOnTop+1)/(Rest+1))) %>%
  filter(log2FC >1.5| log2FC < -1.5)

s1s2_scTransform.markers = as_tibble(read.table(file = paste0(plotDirectory, "scTransformMarkers_snn06.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
s1s2_scTransform.markersSubset <- s1s2_scTransform.markers %>% filter(p_val_adj <0.05, avg_logFC>1)

geneList10X = as_tibble(unique(s1s2_scTransform.markersSubset$gene))
geneListBulk = as_tibble(unique(foldChangeAvgTPM$gene_name))
s1s2_scTransform.markersSubsetBulkSubset = s1s2_scTransform.markersSubset %>% filter(gene %in% geneListBulk$value)

avgTPMselect = foldChangeAvgTPM %>% filter(gene_name %in% geneList10X$value)
upRegulatedGeneList = avgTPMselect %>% filter(log2FC>1.5) %>% select(gene_name)
downRegulatedGeneList = avgTPMselect %>% filter(log2FC< -1.5) %>% select(gene_name)

testCoverageA = s1s2_scTransform.markersSubset %>% group_by(cluster) %>% summarise(fraction = (sum(gene %in% upRegulatedGeneList$gene_name))/(length(gene))) %>% mutate(normalizedFraction = fraction/(max(fraction)))
testCoverage2A = s1s2_scTransform.markersSubset %>% group_by(cluster) %>% summarise(fraction = (sum(gene %in% downRegulatedGeneList$gene_name))/(length(gene))) %>% mutate(normalizedFraction = fraction/(max(fraction)))

###Comment: only 4 genes upregulated compared to the rest. Therefore no similarity analysis possible.


######
# countsTest2 = counts %>% select(Original,colonyTypeFine, isResistant) %>% unique() 
# # #########
# countsPCA = counts %>% filter(!colonyTypeFine == "naive") %>% select(Original, gene_name, tpm)
# countsPCAFilter <- aggregate(countsPCA['tpm'], by=countsPCA[c('Original', 'gene_name')], sum)
# 
# allCounts_rawCountsMatrixNew = countsPCAFilter %>% spread(Original, tpm)
# row.names(allCounts_rawCountsMatrixNew) = allCounts_rawCountsMatrixNew$gene_name
# allCounts_rawCountsMatrixNew$gene_name = NULL
# allCounts_rawCountsMatrixNew[is.na(allCounts_rawCountsMatrixNew)] <- 0
# 
# allCounts_rawCountsMatrixNewFilter = allCounts_rawCountsMatrixNew[rowSums(allCounts_rawCountsMatrixNew[,])> 10,]
# 
# pca_data=prcomp(t(allCounts_rawCountsMatrixNewFilter))
# pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
# 
# condition = colnames(allCounts_rawCountsMatrixNewFilter)
# df_pca_data=data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], sample = colnames(allCounts_rawCountsMatrixNewFilter), condition=condition)
# 
# ggplot(df_pca_data, aes(PC1,PC2, color = sample))+
#   geom_point(size=8)+
#   labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))