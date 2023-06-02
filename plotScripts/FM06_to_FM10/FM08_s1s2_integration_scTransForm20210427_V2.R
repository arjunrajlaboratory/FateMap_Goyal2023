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
dataDirectory <-("~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM08/jointAnalysis/seurat/")
plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM08/jointAnalysis/seuratOutput/"
  
# Load the sample11umPlx dataset
sample11umPlx.data <- Read10X(data.dir = "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM08/Cellranger/run1_sample3/filtered_feature_bc_matrix")
sample21umPlx.data <- Read10X(data.dir = "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM08/Cellranger/run1_sample4/filtered_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
sample11umPlx <- CreateSeuratObject(counts = sample11umPlx.data, project = "10X_Sample1_1uMPLX", min.cells = 3, min.features = 200)
sample21umPlx <- CreateSeuratObject(counts = sample21umPlx.data, project = "10X_Sample2_1uMPLX", min.cells = 3, min.features = 200)

sample11umPlx[["percent.mt"]] <- PercentageFeatureSet(object = sample11umPlx, pattern = "^MT-")
sample21umPlx[["percent.mt"]] <- PercentageFeatureSet(object = sample21umPlx, pattern = "^MT-")

#####
VlnPlot(object = sample11umPlx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = sample11umPlx, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(object = sample11umPlx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
VlnPlot(object = sample21umPlx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot3 <- FeatureScatter(object = sample21umPlx, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot4 <- FeatureScatter(object = sample21umPlx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot3,plot4))
CombinePlots(plots = list(plot3,plot1))

sample11umPlx <- subset(x = sample11umPlx, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 18)
sample21umPlx <- subset(x = sample21umPlx, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 18)

s1s2_scTransform <- merge(sample11umPlx, y = sample21umPlx, add.cell.ids = c("S1", "S2"), project = "S1S2")


s1s2_scTransform.list <- SplitObject(s1s2_scTransform, split.by = "orig.ident")
s1s2_scTransform.list <- s1s2_scTransform.list[c("10X_Sample1_1uMPLX", "10X_Sample2_1uMPLX")]

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
s1s2_scTransform.integrated <- RunUMAP(s1s2_scTransform.integrated, dims = 1:50)

saveRDS(s1s2_scTransform.integrated, file = paste0(dataDirectory,'s1s2ScTransform_50pcs_filter.rds'))


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