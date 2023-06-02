#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2020_FM02/jointAnalysis/seurat/hg19/s1s3s4_Scanoramahg.R
##############################################

library(ggplot2)
library(Seurat)
library(Matrix)
library(stringr)
library(readr)
library(here)
library(fitdistrplus)
library(dplyr)
library(data.table)
options(future.globals.maxSize = 10000 * 1024^2)

#####restart R to run this command first
library(reticulate)
path_to_python <- "/Users/yogesh/opt/anaconda3/"
use_python(path_to_python)
reticulate::use_condaenv(condaenv = "Scanorama")
scanorama = reticulate::import(module = "scanorama")
#####restart R to run this command first ^^

dataDirectory <-("/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM02/")
home1Directory <-("/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM02/")

# Load the sample11umPlx dataset
sample11uMPlx.data <- Read10X(data.dir = "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/rawData/10X_cellRangerOuts/FM02/1-1uMPLX/filtered_feature_bc_matrix/")
sample35nMT.data <- Read10X(data.dir = "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/rawData/10X_cellRangerOuts/FM02/3-5nMT/filtered_feature_bc_matrix/")
sample45nMT100nMPLX.data <- Read10X(data.dir = "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/rawData/10X_cellRangerOuts/FM02/4-5nMT100nMP/filtered_feature_bc_matrix/")

sample11uMPlx <- CreateSeuratObject(counts = sample11uMPlx.data, project = "Sample1_1uMPLX", min.cells = 3, min.features = 200)
sample35nMT <- CreateSeuratObject(counts = sample35nMT.data, project = "Sample3_5nMT", min.cells = 3, min.features = 200)
sample45nMT100nMPLX <- CreateSeuratObject(counts = sample45nMT100nMPLX.data, project = "Sample4_5nMT1uMPLX", min.cells = 3, min.features = 200)

sample11uMPlx[["percent.mt"]] <- PercentageFeatureSet(object = sample11uMPlx, pattern = "^MT-")
sample35nMT[["percent.mt"]] <- PercentageFeatureSet(object = sample35nMT, pattern = "^MT-")
sample45nMT100nMPLX[["percent.mt"]] <- PercentageFeatureSet(object = sample45nMT100nMPLX, pattern = "^MT-")

##### Use this to find cutoffs for percetMT
VlnPlot(object = sample45nMT100nMPLX, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = sample45nMT100nMPLX, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(object = sample45nMT100nMPLX, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

sample11uMPlx <- subset(x = sample11uMPlx, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 30)
sample35nMT <- subset(x = sample35nMT, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 30)
sample45nMT100nMPLX <- subset(x = sample45nMT100nMPLX, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 30)

s1s3s4_scanoramahg19 <- merge(sample11uMPlx, y = c(sample35nMT,sample45nMT100nMPLX), add.cell.ids = c("S1", "S3", "S4"))


extractDataScanorama <- function(seurat.object, assay = "RNA", slot = "counts", groupingVar = "orig.ident", group_name){
  return(t(as.matrix(GetAssayData(seurat.object, assay = assay, slot = slot)))[colnames(seurat.object)[seurat.object@meta.data[,groupingVar] == group_name],])
}

samples = c("Sample1_1uMPLX","Sample3_5nMT", "Sample4_5nMT1uMPLX")

datasets = list()
gene_list = list()

for (i in 1:length(samples)){
  datasets[[i]] <- extractDataScanorama(s1s3s4_scanoramahg19, group_name = samples[[i]])
  gene_list[[i]] <- rownames(s1s3s4_scanoramahg19)
}
rm(sample11uMPlx)
rm(sample11uMPlx.data)
rm(sample35nMT)
rm(sample35nMT.data)
rm(sample45nMT100nMPLX)
rm(sample45nMT100nMPLX.data)

integrated.corrected.data = scanorama$correct(datasets, gene_list, return_dimred=TRUE, return_dense=TRUE, ds_names = samples, verbose = TRUE)

corrected_scanorama <- t(do.call(rbind, integrated.corrected.data[[2]]))
colnames(corrected_scanorama) <- colnames(s1s3s4_scanoramahg19)
rownames(corrected_scanorama) <- integrated.corrected.data[[3]]
corrected_scanorama_pca <- t(do.call(rbind, integrated.corrected.data[[1]]))
colnames(corrected_scanorama_pca) <- colnames(s1s3s4_scanoramahg19)

# add in assay and format as  seurat object
scanorama_assay <- CreateAssayObject(data = corrected_scanorama)
s1s3s4_scanoramahg19[["scanorama"]] <- scanorama_assay
DefaultAssay(s1s3s4_scanoramahg19) <- "scanorama"

#s1s2_scTransform[['scanorama']]@counts = s1s2_scTransform[['scanorama']]@data
#a = s1s2_scTransform[['scanorama']]@counts
#b = length(is.na(a))
# Preprocess scanorama values and perform PCA
s1s3s4_scanoramahg19 <- FindVariableFeatures(s1s3s4_scanoramahg19, assay = "scanorama", selection.method = "vst", nfeatures = 7000)
all.genes <- rownames(s1s3s4_scanoramahg19)
s1s3s4_scanoramahg19 <- ScaleData(s1s3s4_scanoramahg19, features = all.genes)
s1s3s4_scanoramahg19 <- RunPCA(object = s1s3s4_scanoramahg19, assay = "scanorama", reduction.name = "pca_scanorama")
s1s3s4_scanoramahg19 <- FindNeighbors(object=s1s3s4_scanoramahg19, dims=1:50, reduction = "pca_scanorama", force.recalc = TRUE, graph.name = "scanorama_snn")
s1s3s4_scanoramahg19 <- FindClusters(object=s1s3s4_scanoramahg19,graph.name = "scanorama_snn", resolution=0.8)
s1s3s4_scanoramahg19 <- RunUMAP(object = s1s3s4_scanoramahg19, reduction = "pca_scanorama", dims = 1:50, reduction.name = "umap_scanorama")
clusterUMAP = DimPlot(s1s3s4_scanoramahg19, reduction = "umap_scanorama", label = TRUE, group.by = "scanorama_snn_res.0.8")
sampleUMAP = DimPlot(s1s3s4_scanoramahg19, reduction = "umap_scanorama", group.by = "orig.ident")
sampleSplitUMAP = DimPlot(s1s3s4_scanoramahg19, reduction = "umap_scanorama",slot="RNA.data", split.by = "orig.ident")

saveRDS(s1s3s4_scanoramahg19, file = paste0(dataDirectory,'s1s3s4_scanorama_50pcs_filter_hg19.rds'))

#############################################################################
#DON'T NEED TO RUN BLOCK THIS AGAIN
#############################################################################

s1s3s4Scanorama <- readRDS(paste0(home1Directory,"s1s3s4_scanorama_50pcs_filter_hg19.rds"))

logNormalizedCounts = s1s3s4Scanorama@assays$scanorama@scale.data
logNormalizedCountsRound = round(logNormalizedCounts,4)

#s1s2Scanorama[["scanorama"]]@scale.data

cells_count = colnames(logNormalizedCountsRound) #CellIds with Sample number as prefix
cells_count_cellID = sub("S\\d_", "", cells_count)
cells_count_sampleNum = gsub("[^S134]", "", cells_count)
logNormalizedCountsRound = as_tibble(as.data.frame((t(as.matrix(logNormalizedCountsRound)))))
logNormalizedCountsRound = logNormalizedCountsRound %>% mutate(cellID = cells_count_cellID,
                                                               sampleNum = cells_count_sampleNum)

umapCoordinates = (s1s3s4Scanorama[['umap_scanorama']])@cell.embeddings #UMAP coordinates
cells_UMAP = rownames(umapCoordinates) #CellIds with Sample number as prefix
cells_UMAP_cellID = sub("S\\d_", "", cells_UMAP)
cells_UMAP_sampleNum = gsub("[^S134]", "", cells_UMAP)
umapCoordinates = as_tibble(umapCoordinates) #UMAP coordinates in tibble
umapCoordinates = umapCoordinates %>% mutate(cellID = cells_UMAP_cellID,
                                             sampleNum = cells_UMAP_sampleNum)

s1s3s4Scanorama <- FindClusters(object=s1s3s4Scanorama,graph.name = "scanorama_snn", resolution=0.5)

umapClusters = (s1s3s4Scanorama[['scanorama_snn_res.0.5']])
cells_Clusters = rownames(umapClusters) #CellIds with Sample number as prefix
cells_Clusters_cellID = sub("S\\d_", "", cells_Clusters)
cells_Clusters_sampleNum = gsub("[^S134]", "", cells_Clusters)
umapClusters = as_tibble(umapClusters)
umapClusters = umapClusters %>% mutate(cellID = cells_Clusters_cellID,
                                       sampleNum = cells_Clusters_sampleNum)

write.table(umapCoordinates, file=paste0(home1Directory,'umapCoordinates_s1s3s4Scanorama_50pcs_filter_hg19.tsv'), col.names = TRUE, sep='\t')
#write.table(logNormalizedCounts, file=paste0(home1Directory,'logNormalizedCounts_s1s2Scanorama_50pcs_filter.tsv'), col.names = TRUE, sep='\t')
write.table(umapClusters, file=paste0(home1Directory,'umapClusters_s1s3s4Scanorama_snn05_50pcs_filter_hg19.tsv'), col.names = TRUE, sep='\t')
write.table(logNormalizedCountsRound, file=paste0(home1Directory,'logNormalizedCounts_s1s3s4Scanorama_50pcs_filterRound_hg19.tsv'), col.names = TRUE, sep='\t')

rm(s1s3s4Scanorama)
rm(umapCoordinates)
rm(umapClusters)
rm(logNormalizedCounts)
rm(logNormalizedCountsRound)

#############################################################################
#DON'T NEED TO RUN BLOCK THIS AGAIN
#############################################################################