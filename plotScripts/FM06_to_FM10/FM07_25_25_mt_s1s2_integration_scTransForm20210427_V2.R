#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2019_FM01/SubmissionScripts/FINAL_COPIED/finalFiguresPaper/s1s2_integration_scTransForm20210427.R
##############################################
########################################################################################
#Description: Processing and analysis of FM07
#Author: Yogesh Goyal
#Edited by: Maalavika Pillai
#Version: 2.1
#Edited on: 10/11/22
########################################################################################

library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggrepel)
options(future.globals.maxSize = 4000 * 1024^2)

#########
# workflow adapted from https://satijalab.org/seurat/v3.0/integration.html
#########
 rm(sample11umPlx)
rm(sample11umPlx.data)
rm(sample21umPlx)
rm(sample21umPlx.data)
rm(s1s2_scTransform)
rm(s1s2_scTransform.list)
rm(s1s2_scTransform.anchors)

#plots <- DimPlot(s1s2_scTransform.integrated, reduction = "umap", group.by = "orig.ident", combine = FALSE)



###Functionalize Count normalized matrix
logNormalizedCounts = s1s2_scTransform.integrated[['SCT']]@data #normalized log counts matrix (for counts only, replace "data" by "counts")
normalizedCounts = s1s2_scTransform.integrated[['SCT']]@counts #normalized counts matrix
cells_count = colnames(normalizedCounts) #CellIds with Sample number as prefix
cells_count_cellID = sub("S\\d_", "", cells_count)
cells_count_sampleNum = gsub("[^S12]", "", cells_count)
logNormalizedCounts = as_tibble(as.data.frame((t(as.matrix(logNormalizedCounts)))))
logNormalizedCounts = logNormalizedCounts %>% mutate(cellID = cells_count_cellID,
                                                     sampleNum = cells_count_sampleNum)
normalizedCounts = as_tibble(as.data.frame((t(as.matrix(normalizedCounts)))))
normalizedCounts = normalizedCounts %>% mutate(cellID = cells_count_cellID,
                                                     sampleNum = cells_count_sampleNum)

write.table(logNormalizedCounts, file=paste0(dataDirectory,'logNormalizedSCTCountsS1S2_50pcs_filter.tsv'), col.names = TRUE, sep='\t')
write.table(normalizedCounts, file=paste0(dataDirectory,'normalizedSCTCountsS1S2_50pcs_filter.tsv'), col.names = TRUE, sep='\t')

#test =  (s1s2_scTransform.integrated[['orig.ident']])
umapCoordinates = (s1s2_scTransform.integrated[['umap']])@cell.embeddings #UMAP coordinates
cells_UMAP = rownames(umapCoordinates) #CellIds with Sample number as prefix
cells_UMAP_cellID = sub("S\\d_", "", cells_UMAP)
cells_UMAP_sampleNum = gsub("[^S12]", "", cells_UMAP)
umapCoordinates = as_tibble(umapCoordinates) #UMAP coordinates in tibble
umapCoordinates = umapCoordinates %>% mutate(cellID = cells_UMAP_cellID,
                                             sampleNum = cells_UMAP_sampleNum)
write.table(umapCoordinates, file=paste0(dataDirectory,'umapCoordinatesSCT_S1S2_50pcs_filter.tsv'), col.names = TRUE, sep='\t')