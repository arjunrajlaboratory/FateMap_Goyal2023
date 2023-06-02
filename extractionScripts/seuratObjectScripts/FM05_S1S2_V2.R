#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2021_FM05_06/jointAnalysis/seurat/FM05_S1S2.R
##############################################

library(Seurat)
library(ggplot2)
library(tidyverse)
options(future.globals.maxSize = 4000 * 1024^2)

#########
# workflow adapted from https://satijalab.org/seurat/v3.0/integration.html
#########

dataDirectory <-("/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM05/")
home1Directory <-("/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM05/")
# Load the sample11umPlx dataset
sample1250nM.data <- Read10X(data.dir = "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/rawData/10X_cellRangerOuts/FM05/FM05-250nMPLX-A1/filtered_feature_bc_matrix/")
sample2250nM.data <- Read10X(data.dir = "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/rawData/10X_cellRangerOuts/FM05/FM05-250nMPLX-A2/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
sample1250nM <- CreateSeuratObject(counts = sample1250nM.data, project = "10X_Sample1_250nM", min.cells = 3, min.features = 200)
sample2250nM <- CreateSeuratObject(counts = sample2250nM.data, project = "10X_Sample2_250nM", min.cells = 3, min.features = 200)

sample1250nM[["percent.mt"]] <- PercentageFeatureSet(object = sample1250nM, pattern = "^MT-")
sample2250nM[["percent.mt"]] <- PercentageFeatureSet(object = sample2250nM, pattern = "^MT-")


#####
VlnPlot(object = sample2250nM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = sample2250nM, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(object = sample2250nM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))

sample1250nM <- subset(x = sample1250nM, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)
sample2250nM <- subset(x = sample2250nM, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)

s1s2_scTransform <- merge(sample1250nM, y = sample2250nM, add.cell.ids = c("S1", "S2"), project = "S1S2")
s1s2_scTransform

s1s2_scTransform.list <- SplitObject(s1s2_scTransform, split.by = "orig.ident")
s1s2_scTransform.list <- s1s2_scTransform.list[c("10X_Sample1_250nM", "10X_Sample2_250nM")]

for (i in 1:length(s1s2_scTransform.list)) {
  s1s2_scTransform.list[[i]] <- SCTransform(s1s2_scTransform.list[[i]], verbose = FALSE)
}

s1s2_scTransform.features <- SelectIntegrationFeatures(object.list = s1s2_scTransform.list, nfeatures = 7000)
s1s2_scTransform.list <- PrepSCTIntegration(object.list = s1s2_scTransform.list, anchor.features = s1s2_scTransform.features, 
                                            verbose = FALSE)

s1s2_scTransform.anchors <- FindIntegrationAnchors(object.list = s1s2_scTransform.list, normalization.method = "SCT", 
                                                   anchor.features = s1s2_scTransform.features, verbose = FALSE)
s1s2_scTransform.integrated <- IntegrateData(anchorset = s1s2_scTransform.anchors, normalization.method = "SCT", 
                                             verbose = FALSE)

s1s2_scTransform.integrated <- RunPCA(s1s2_scTransform.integrated, verbose = FALSE)
s1s2_scTransform.integrated <- FindNeighbors(object=s1s2_scTransform.integrated, dims=1:50, verbose = FALSE)
s1s2_scTransform.integrated <- FindClusters(object=s1s2_scTransform.integrated, resolution = 0.5, verbose = FALSE)
s1s2_scTransform.integrated <- RunUMAP(s1s2_scTransform.integrated, dims = 1:50)

saveRDS(s1s2_scTransform.integrated, file = paste0(dataDirectory,'s1s2250nMWM983B_ScTransform_50pcs_filter.rds'))

s1s2ScTransform <- readRDS(paste0(home1Directory,"s1s2250nMWM983B_ScTransform_50pcs_filter.rds"))

logNormalizedCounts = s1s2ScTransform[['SCT']]@data #normalized log counts matrix (for counts only, replace "data" by "counts")
normalizedCounts = s1s2ScTransform[['SCT']]@counts #normalized counts matrix
cells_count = colnames(normalizedCounts) #CellIds with Sample number as prefix
cells_count_cellID = sub("S\\d_", "", cells_count)
cells_count_sampleNum = gsub("[^S12]", "", cells_count)
logNormalizedCounts = as_tibble(as.data.frame((t(as.matrix(logNormalizedCounts)))))
logNormalizedCounts = logNormalizedCounts %>% mutate(cellID = cells_count_cellID,
                                                     sampleNum = cells_count_sampleNum)
normalizedCounts = as_tibble(as.data.frame((t(as.matrix(normalizedCounts)))))
normalizedCounts = normalizedCounts %>% mutate(cellID = cells_count_cellID,
                                               sampleNum = cells_count_sampleNum)

umapCoordinates = (s1s2ScTransform[['umap']])@cell.embeddings #UMAP coordinates
cells_UMAP = rownames(umapCoordinates) #CellIds with Sample number as prefix
cells_UMAP_cellID = sub("S\\d_", "", cells_UMAP)
cells_UMAP_sampleNum = gsub("[^S12]", "", cells_UMAP)
umapCoordinates = as_tibble(umapCoordinates) #UMAP coordinates in tibble
umapCoordinates = umapCoordinates %>% mutate(cellID = cells_UMAP_cellID,
                                             sampleNum = cells_UMAP_sampleNum)

write.table(logNormalizedCounts, file=paste0(home1Directory,'logNormalizedSCTCountsS1S2_FM05_50pcs_filter.tsv'), col.names = TRUE, sep='\t')
write.table(normalizedCounts, file=paste0(home1Directory,'normalizedSCTCountsS1S2_FM05_50pcs_filter.tsv'), col.names = TRUE, sep='\t')
write.table(umapCoordinates, file=paste0(home1Directory,'umapCoordinatesSCT_S1S2_FM05_50pcs_filter.tsv'), col.names = TRUE, sep='\t')