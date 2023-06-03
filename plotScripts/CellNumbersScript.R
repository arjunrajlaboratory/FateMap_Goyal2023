library(ggplot2)
library(tidyverse)
library(reshape2)
library(svglite)
library(Seurat)

#### Memory across sampels_using data from s1s2_integration_scTransForm
home1Directory <-"/Users/ygy1258/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM01/"
umapCoordinates = as_tibble(read.table(file = paste0(home1Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum == "S1")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum == "S2")


home2Directory <-"/Users/ygy1258/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM02/"
umapCoordinates = as_tibble(read.table(file = paste0(home2Directory, "umapCoordinates_s1s2Scanorama_50pcs_filter_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum == "S1")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum == "S2")

umapCoordinates = as_tibble(read.table(file = paste0(home2Directory, "umapCoordinates_s1s3s4Scanorama_50pcs_filter_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum == "S1")
umapCoordinatesS3 = umapCoordinates %>% filter(sampleNum == "S3")
umapCoordinatesS4 = umapCoordinates %>% filter(sampleNum == "S4")

home3Directory <-"/Users/ygy1258/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM03/"
umapCoordinates = as_tibble(read.table(file = paste0(home3Directory, "umapCoordinates_s1s2Scanorama_50pcs_filter_hg19.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum == "S1")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum == "S2")


home4Directory <-"/Users/ygy1258/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM04/"
umapCoordinates = as_tibble(read.table(file = paste0(home4Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum == "S1")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum == "S2")

home5Directory <-"/Users/ygy1258/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM05/"
umapCoordinates = as_tibble(read.table(file = paste0(home5Directory, "umapCoordinatesSCT_S1S2_FM05_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum == "S1")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum == "S2")

home5Directory <-"/Users/ygy1258/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM05/"
umapCoordinates = as_tibble(read.table(file = paste0(home5Directory, "umapCoordinatesSCT_S1S2_FM05_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum == "S1")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum == "S2")

##check
home6Directory <-"/Users/ygy1258/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2021_FM06/extracted/FM06/"
fm06 <- readRDS(paste0(home6Directory,"s1s2Naive_ScTransform_50pcs_filter.rds"))
fm06.updated = UpdateSeuratObject(object = fm06)
umapCoordinates = (fm06.updated[['umap']])@cell.embeddings #UMAP coordinates
cells_UMAP = rownames(umapCoordinates) #CellIds with Sample number as prefix
cells_UMAP_cellID = sub("S\\d_", "", cells_UMAP)
cells_UMAP_sampleNum = gsub("[^S12]", "", cells_UMAP)
umapCoordinates = as_tibble(umapCoordinates) #UMAP coordinates in tibble
umapCoordinates = umapCoordinates %>% mutate(cellID = cells_UMAP_cellID,
                                             sampleNum = cells_UMAP_sampleNum)
write.table(umapCoordinates, file=paste0(home6Directory,'umapCoordinatesSCT_S1S2_50pcs_filter.tsv'), col.names = TRUE, sep='\t')
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum == "S1")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum == "S2")

home7Directory <-"/Users/ygy1258/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM07/jointAnalysis/seurat/"
umapCoordinates = as_tibble(read.table(file = paste0(home7Directory, "25_15_mt_umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum == "S11")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum == "S21")

home8Directory <-"/Users/ygy1258/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM08/jointAnalysis/seurat/"
umapCoordinates = as_tibble(read.table(file = paste0(home8Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum == "S11")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum == "S21")

home9Directory <-"/Users/ygy1258/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM09/jointAnalysis/seurat/"
umapCoordinates = as_tibble(read.table(file = paste0(home9Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum == "S11")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum == "S21")

home10Directory <-"/Users/ygy1258/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM10/jointAnalysis/seurat/"
umapCoordinates = as_tibble(read.table(file = paste0(home10Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum == "S11")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum == "S21")