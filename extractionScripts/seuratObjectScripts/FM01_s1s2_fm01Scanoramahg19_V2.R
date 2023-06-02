#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2019_FM01/SubmissionScripts/FINAL_COPIED/finalFiguresPaper/s1s2_fm01Scanoramahg19.R
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

dataDirectory = "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM01/"
plotDirectory <- "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/plotData/analysisPlots/"

# Load the sample11umPlx dataset
sample11umPlx.data <- Read10X(data.dir = "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/rawData/10X_cellRangerOuts/FM01/1_1uMPLX/filtered_feature_bc_matrix/1_1uMPLX/")
sample21umPlx.data <- Read10X(data.dir = "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/rawData/10X_cellRangerOuts/FM01/2_1uMPLX/filtered_feature_bc_matrix/2_1uMPLX/")

# Initialize the Seurat object with the raw (non-normalized data).
sample11umPlx <- CreateSeuratObject(counts = sample11umPlx.data, project = "10X_Sample1_1uMPLX", min.cells = 3, min.features = 200)
sample21umPlx <- CreateSeuratObject(counts = sample21umPlx.data, project = "10X_Sample2_1uMPLX", min.cells = 3, min.features = 200)

sample11umPlx[["percent.mt"]] <- PercentageFeatureSet(object = sample11umPlx, pattern = "^MT-")
sample21umPlx[["percent.mt"]] <- PercentageFeatureSet(object = sample21umPlx, pattern = "^MT-")

sample11umPlx <- subset(x = sample11umPlx, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)
sample21umPlx <- subset(x = sample21umPlx, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)

s1s2_scTransform <- merge(sample11umPlx, y = sample21umPlx, add.cell.ids = c("S1", "S2"))


extractDataScanorama <- function(seurat.object, assay = "RNA", slot = "counts", groupingVar = "orig.ident", group_name){
  return(t(as.matrix(GetAssayData(seurat.object, assay = assay, slot = slot)))[colnames(seurat.object)[seurat.object@meta.data[,groupingVar] == group_name],])
}

samples = c("10X_Sample1_1uMPLX","10X_Sample2_1uMPLX")

datasets = list()
gene_list = list()

for (i in 1:length(samples)){
  datasets[[i]] <- extractDataScanorama(s1s2_scTransform, group_name = samples[[i]])
  gene_list[[i]] <- rownames(s1s2_scTransform)
}
rm(sample11umPlx)
rm(sample11umPlx.data)
rm(sample21umPlx)
rm(sample21umPlx.data)

integrated.corrected.data = scanorama$correct(datasets, gene_list, return_dimred=TRUE, return_dense=TRUE, ds_names = samples, verbose = TRUE)

corrected_scanorama <- t(do.call(rbind, integrated.corrected.data[[2]]))
colnames(corrected_scanorama) <- colnames(s1s2_scTransform)
rownames(corrected_scanorama) <- integrated.corrected.data[[3]]
corrected_scanorama_pca <- t(do.call(rbind, integrated.corrected.data[[1]]))
colnames(corrected_scanorama_pca) <- colnames(s1s2_scTransform)

# add in assay and format as  seurat object
scanorama_assay <- CreateAssayObject(data = corrected_scanorama)
s1s2_scTransform[["scanorama"]] <- scanorama_assay
DefaultAssay(s1s2_scTransform) <- "scanorama"

s1s2_scTransform <- FindVariableFeatures(s1s2_scTransform, assay = "scanorama", selection.method = "vst", nfeatures = 7000)
all.genes <- rownames(s1s2_scTransform)
s1s2_scTransform <- ScaleData(s1s2_scTransform, features = all.genes)
s1s2_scTransform <- RunPCA(object = s1s2_scTransform, assay = "scanorama", reduction.name = "pca_scanorama")
s1s2_scTransform <- FindNeighbors(object=s1s2_scTransform, dims=1:50, reduction = "pca_scanorama", force.recalc = TRUE, graph.name = "scanorama_snn")
s1s2_scTransform <- FindClusters(object=s1s2_scTransform,graph.name = "scanorama_snn", resolution=0.8)
s1s2_scTransform <- RunUMAP(object = s1s2_scTransform, reduction = "pca_scanorama", dims = 1:50, reduction.name = "umap_scanorama")
clusterUMAP = DimPlot(s1s2_scTransform, reduction = "umap_scanorama", label = TRUE, group.by = "scanorama_snn_res.0.8")
sampleUMAP = DimPlot(s1s2_scTransform, reduction = "umap_scanorama", group.by = "orig.ident")
sampleSplitUMAP = DimPlot(s1s2_scTransform, reduction = "umap_scanorama",slot="RNA.data", split.by = "orig.ident")

saveRDS(s1s2_scTransform, file = paste0(dataDirectory,'s1s2_fm01_scanorama_50pcs_filter_MT26_hg19.rds'))

s1s2Scanorama <- readRDS(paste0(dataDirectory,"s1s2_fm01_scanorama_50pcs_filter_MT26_hg19.rds"))

clusterUMAP_08 = DimPlot(s1s2Scanorama, reduction = "umap_scanorama", label = FALSE, group.by = "seurat_clusters") + NoLegend()
clusterUMAP_08WCluster = DimPlot(s1s2Scanorama, reduction = "umap_scanorama", label = TRUE, label.size = 8, group.by = "seurat_clusters") + NoLegend()
sampleUMAP = DimPlot(s1s2Scanorama, reduction = "umap_scanorama", group.by = "orig.ident") + NoLegend()
sampleSplitUMAP_08 = DimPlot(s1s2Scanorama, reduction = "umap_scanorama",slot="RNA.data", group.by = "seurat_clusters", split.by = "orig.ident", label = FALSE) + NoLegend()

ggsave(clusterUMAP_08, file = paste0(plotDirectory, 'FM01_s1s2_clusterUMAP_onlyScanorma08_hg19.svg'), width = 6, height = 5.321)
ggsave(clusterUMAP_08WCluster, file = paste0(plotDirectory, 'FM01_s1s2_clusterUMAPWNUMBERS_onlyScanorma08_hg19.svg'), width = 6, height = 5.321)
ggsave(sampleUMAP, file = paste0(plotDirectory, 'FM01_s1s2_sampleUMAP_onlyScanorma_hg19.svg'), width = 6, height = 5.321)
ggsave(sampleSplitUMAP_08, file = paste0(plotDirectory, 'FM01_s1s2_sampleSplitUMAP_onlyScanorma08_hg19.svg'), width = 12, height = 5.822)
