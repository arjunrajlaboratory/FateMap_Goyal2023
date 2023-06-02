#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/RNASeq/20200205_RNASeq/Downstream_R_Analysis/TPM_Analysis.R
##############################################

####script to analyze/heatmap for fast vs slow invadig cells based on individual colonies
library(dplyr)
library(ggplot2)
library(gplots)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)

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

#counts table for all XX samples

counts = as_tibble(read.table(file = paste0(dataDirectory,"deRandomized_BulkColonySeq_20200220_meltedData_hg19_normalizedRevised.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))

counts1 = counts %>% dplyr::select(tpm, gene_name,Original)

countsFilter = counts1 %>% filter(Original %in% c("sample_38", "sample_40","sample_42","sample_43"))
countsFilter1 <- aggregate(countsFilter['tpm'], by=countsFilter[c('Original', 'gene_name')], sum) ##to get rid of same genes assigned multiple times from htseq

allCounts_rawCountsMatrixNew = unique(countsFilter1) %>%
  spread(Original, tpm)

allCounts_rawCountsMatrixNew[is.na(allCounts_rawCountsMatrixNew)] <- 0

#Make gene_id the rowname: prevents losing the gene_id during column selection
row.names(allCounts_rawCountsMatrixNew) = allCounts_rawCountsMatrixNew$gene_name
allCounts_rawCountsMatrixNew$gene_name = NULL

#log of the matrix (1 added to avoid log(0) errors). Default log is log_e
logTPM <- log2(allCounts_rawCountsMatrixNew + 1)

#eliminate low tpm values: done by calculating sum across each row. 
logTPM1 <- logTPM[!(rowSums(logTPM)<=3),]    

var_genes <- apply(logTPM1, 1, var) ##variance
head(var_genes)

select_var <- names(sort(var_genes, decreasing=TRUE))[1:200]
head(select_var)

highly_variable_lcpm <-  as.matrix(logTPM1[select_var,])

###log2FC analysis
avgTPM = allCounts_rawCountsMatrixNew %>% rownames_to_column(var = "gene_name") %>% 
  mutate(avgSlow = (sample_38+sample_40)/2, avgFast = (sample_42+sample_43)/2) %>%
  filter((avgFast+avgSlow)/2 > 10) %>% mutate(log2FC = log2((avgFast+1)/(avgSlow+1))) %>%
  filter(log2FC >1.5| log2FC < -1.5)

s1s2_scTransform.markers = as_tibble(read.table(file = paste0(plotDirectory, "scTransformMarkers_snn06.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
s1s2_scTransform.markersSubset <- s1s2_scTransform.markers %>% filter(p_val_adj <0.05, avg_logFC>1)

geneList10X = as_tibble(unique(s1s2_scTransform.markersSubset$gene))
geneListBulk = as_tibble(unique(avgTPM$gene_name))

s1s2_scTransform.markersSubsetBulkSubset = s1s2_scTransform.markersSubset %>% filter(gene %in% geneListBulk$value)
avgTPMselect = avgTPM %>% filter(gene_name %in% geneList10X$value)
upRegulatedGeneList = avgTPMselect %>% filter(log2FC>1.5) %>% select(gene_name)
downRegulatedGeneList = avgTPMselect %>% filter(log2FC< -1.5) %>% select(gene_name)

avgTPMselectFast = avgTPMselect %>% select(gene_name, log2FC) %>% filter(log2FC>0)
write.table(avgTPMselectFast, file=paste0(plotDirectory,'fastMigratingGeneList.tsv'), col.names = TRUE, sep='\t')

avgTPMselectSlow = avgTPMselect %>% select(gene_name, log2FC) %>% filter(log2FC<0) %>% mutate(log2FC = -log2FC)
write.table(avgTPMselectSlow, file=paste0(plotDirectory,'slowMigratingGeneList.tsv'), col.names = TRUE, sep='\t')

testCoverageA = s1s2_scTransform.markersSubset %>% group_by(cluster) %>% summarise(fraction = (sum(gene %in% upRegulatedGeneList$gene_name))/(length(gene))) %>% mutate(normalizedFraction = fraction/(max(fraction)))
testCoverage2A = s1s2_scTransform.markersSubset %>% group_by(cluster) %>% summarise(fraction = (sum(gene %in% downRegulatedGeneList$gene_name))/(length(gene))) %>% mutate(normalizedFraction = fraction/(max(fraction)))

#testCoverage = s1s2_scTransform.markersSubset %>% group_by(cluster) %>% summarise(fraction = (sum(gene %in% upRegulatedGeneList$gene_name))) %>% mutate(normalizedFraction = fraction/(max(fraction)))
#testCoverage2 = s1s2_scTransform.markersSubset %>% group_by(cluster) %>% summarise(fraction = (sum(gene %in% downRegulatedGeneList$gene_name))) %>% mutate(normalizedFraction = fraction/(max(fraction)))


test = s1s2_scTransform.markersSubset %>%  group_by(cluster) %>% summarise(totalGenes = length(gene))

    
umapCoordinates = as_tibble(read.table(file = paste0(home1Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapClusters = as_tibble(read.table(file = paste0(home1Directory, "umapClusters_s1s2Scanorama_50pcs_filter_snn0_6.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t")) %>% rename(cluster = seurat_clusters)

umapClustersFractionUp = inner_join(umapClusters,testCoverageA, by = "cluster" )
umapClustersFractionDown = inner_join(umapClusters,testCoverage2A, by = "cluster" )

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
ggsave(plot, file = paste0(plotDirectory, 'phenotypeFastMigratingMap.svg'), width = 6, height = 6.429)

plot <- ggplot(umapCoordinates, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes_string(color = umapClustersFractionDown$normalizedFraction), size = 1, shape = 16) +
  scale_color_gradient(low = "gray93", high = "red4") +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color =  "Scaled\n Log UMI counts")
ggsave(plot, file = paste0(plotDirectory, 'phenotypeSlowMigratingMap.svg'), width = 6, height = 6.429)

####HeatMap from Variable Genes
avgTPMselectHeatMap = avgTPMselect %>% select(-avgSlow, -avgFast, -log2FC)
allTableSpread = as.matrix(avgTPMselectHeatMap[,2:5])
rownames(allTableSpread) = avgTPMselectHeatMap$gene_name
allTableSpread <- t(allTableSpread)

breaks = seq(-2,2, by = 0.05) #hard code scale from -2 to 2
col <- colorRampPalette(brewer.pal(n=10, "RdYlBu"))(length(breaks)) #set colors
col <- rev(col)
#plot
p1 <- pheatmap(allTableSpread, cluster_rows = T, cluster_cols = T,
               angle_col = 90, fontsize_col = 4,fontsize_row = 4, scale = 'column',
               color = col)
ggsave(p1, file = paste0(plotDirectory, 'bulk_FastSlow_Heatmap.svg'), width = 14, height = 4)

###############All Trametinib samples together
countsFilter = counts1 %>% filter(Original %in% c("sample_38", "sample_40","sample_42","sample_43", "sample_50", "sample_48", "sample_49", "sample_41"))
countsFilter1 <- aggregate(countsFilter['tpm'], by=countsFilter[c('Original', 'gene_name')], sum)

allCounts_rawCountsMatrixNew = unique(countsFilter1) %>%
  spread(Original, tpm)

allCounts_rawCountsMatrixNew[is.na(allCounts_rawCountsMatrixNew)] <- 0

#Make gene_id the rowname: prevents losing the gene_id during column selection
row.names(allCounts_rawCountsMatrixNew) = allCounts_rawCountsMatrixNew$gene_name
allCounts_rawCountsMatrixNew$gene_name = NULL

#log of the matrix (1 added to avoid log(0) errors). Default log is log_e
logTPM <- log2(allCounts_rawCountsMatrixNew + 1)

#eliminate low tpm values: done by calculating sum across each row. 
logTPM1 <- logTPM[!(rowSums(logTPM)<=3),]    

var_genes <- apply(logTPM1, 1, var) ##variance
head(var_genes)

select_var <- names(sort(var_genes, decreasing=TRUE))[1:200]
head(select_var)

highly_variable_lcpm <-  as.matrix(logTPM1[select_var,])

#create settings for color scaling  
breaks = seq(-2,2, by = 0.05) #hard code scale from -2 to 2
col <- colorRampPalette(brewer.pal(n=10, "RdYlBu"))(length(breaks)) #set colors
col <- rev(col)
#plot
p1 <- pheatmap(highly_variable_lcpm, cluster_rows = T, cluster_cols = T,
               angle_col = 45, fontsize_col = 7,fontsize_row = 5, scale = 'row',
               color = col)
ggsave(p1, file = paste0(plotDirectory, 'bulk_FastSlowMedium_Heatmap.svg'), width = 4, height = 9)

################
#create settings for color scaling  OLD heatmap
# breaks = seq(-2,2, by = 0.05) #hard code scale from -2 to 2
# col <- colorRampPalette(brewer.pal(n=10, "RdYlBu"))(length(breaks)) #set colors
# col <- rev(col)
# #plot
# p1 <- pheatmap(highly_variable_lcpm, cluster_rows = T, cluster_cols = T,
#                angle_col = 45, fontsize_col = 7,fontsize_row = 7, scale = 'row',
#                color = col)
# ggsave(p1, file = paste0(plotDirectory, 'bulk_FastSlow_Heatmap200Var.svg'), width = 4, height = 7)