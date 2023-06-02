#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2020_FM03/jointAnalysis/seurat/hg19/s1s2_fm03_scanorama_V1.R
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

library(reticulate)
path_to_python <- "/Users/yogesh/opt/anaconda3/"
use_python(path_to_python)
reticulate::use_condaenv(condaenv = "Scanorama")
scanorama = reticulate::import(module = "scanorama")

dataDirectory <-("/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM03/")
home1Directory <-("/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM03/")
#plotDirectory <- ("/Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2020_FM03/jointAnalysis/seuratOutputData/plots/")
#####################################################
###NOTE: sample names are named as sample11uMPlx and sample2100nMPlx BUT they are, in reality, DMSO, and DOT1L
#####################################################
# Load the sample11umPlx dataset
sample11uMPlx.data <- Read10X(data.dir = "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/rawData/10X_cellRangerOuts/FM03/DMSO-FM3-1uMPLX/filtered_feature_bc_matrix/") #DMSO
sample2100nMPlx.data <- Read10X(data.dir = "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/rawData/10X_cellRangerOuts/FM03/DOT1Li-FM3-1uMPLX/filtered_feature_bc_matrix/") #DOT1L

sample11uMPlx <- CreateSeuratObject(counts = sample11uMPlx.data, project = "Sample1", min.cells = 3, min.features = 200)
sample2100nMPlx <- CreateSeuratObject(counts = sample2100nMPlx.data, project = "Sample2", min.cells = 3, min.features = 200)

sample11uMPlx[["percent.mt"]] <- PercentageFeatureSet(object = sample11uMPlx, pattern = "^MT-")
sample2100nMPlx[["percent.mt"]] <- PercentageFeatureSet(object = sample2100nMPlx, pattern = "^MT-")

#####
#VlnPlot(object = sample11uMPlx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(object = sample2100nMPlx, feature1 = "nCount_RNA", feature2 = "percent.mt") 
#plot2 <- FeatureScatter(object = sample2100nMPlx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1,plot2))
sample11uMPlx <- subset(x = sample11uMPlx, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 30) ##checked
sample2100nMPlx <- subset(x = sample2100nMPlx, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 30) ##checked

s1s2_scTransform <- merge(sample11uMPlx, y = sample2100nMPlx, add.cell.ids = c("S1", "S2"))


extractDataScanorama <- function(seurat.object, assay = "RNA", slot = "counts", groupingVar = "orig.ident", group_name){
  return(t(as.matrix(GetAssayData(seurat.object, assay = assay, slot = slot)))[colnames(seurat.object)[seurat.object@meta.data[,groupingVar] == group_name],])
}

samples = c("Sample1","Sample2")

datasets = list()
gene_list = list()

for (i in 1:length(samples)){
  datasets[[i]] <- extractDataScanorama(s1s2_scTransform, group_name = samples[[i]])
  gene_list[[i]] <- rownames(s1s2_scTransform)
}
rm(sample11uMPlx)
rm(sample11uMPlx.data)
rm(sample2100nMPlx)
rm(sample2100nMPlx.data)

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

#s1s2_scTransform[['scanorama']]@counts = s1s2_scTransform[['scanorama']]@data
#a = s1s2_scTransform[['scanorama']]@counts
#b = length(is.na(a))
# Preprocess scanorama values and perform PCA
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

saveRDS(s1s2_scTransform, file = paste0(dataDirectory,'s1s2_fm03_scanorama_50pcs_filter_MT30_hg19.rds'))


rm(s1s2_scTransform)
#############################################################################
#DON'T NEED TO RUN BLOCK THIS AGAIN
#############################################################################

s1s2Scanorama <- readRDS(paste0(dataDirectory,"s1s2_fm03_scanorama_50pcs_filter_MT30_hg19.rds"))
s1s2Scanorama <- FindClusters(object=s1s2Scanorama,graph.name = "scanorama_snn", resolution=0.8)

logNormalizedCounts = s1s2Scanorama@assays$scanorama@scale.data
logNormalizedCountsRound = round(logNormalizedCounts,4)

#s1s2Scanorama[["scanorama"]]@scale.data

cells_count = colnames(logNormalizedCountsRound) #CellIds with Sample number as prefix
cells_count_cellID = sub("S\\d_", "", cells_count)
cells_count_sampleNum = gsub("[^S12]", "", cells_count)
logNormalizedCountsRound = as_tibble(as.data.frame((t(as.matrix(logNormalizedCountsRound)))))
logNormalizedCountsRound = logNormalizedCountsRound %>% mutate(cellID = cells_count_cellID,
                                                               sampleNum = cells_count_sampleNum)

umapCoordinates = (s1s2Scanorama[['umap_scanorama']])@cell.embeddings #UMAP coordinates
cells_UMAP = rownames(umapCoordinates) #CellIds with Sample number as prefix
cells_UMAP_cellID = sub("S\\d_", "", cells_UMAP)
cells_UMAP_sampleNum = gsub("[^S12]", "", cells_UMAP)
umapCoordinates = as_tibble(umapCoordinates) #UMAP coordinates in tibble
umapCoordinates = umapCoordinates %>% mutate(cellID = cells_UMAP_cellID,
                                             sampleNum = cells_UMAP_sampleNum)

umapClusters = (s1s2Scanorama[['seurat_clusters']])
cells_Clusters = rownames(umapClusters) #CellIds with Sample number as prefix
cells_Clusters_cellID = sub("S\\d_", "", cells_Clusters)
cells_Clusters_sampleNum = gsub("[^S12]", "", cells_Clusters)
umapClusters = as_tibble(umapClusters)
umapClusters = umapClusters %>% mutate(cellID = cells_Clusters_cellID,
                                       sampleNum = cells_Clusters_sampleNum)

write.table(umapCoordinates, file=paste0(home1Directory,'umapCoordinates_s1s2Scanorama_50pcs_filter_hg19.tsv'), col.names = TRUE, sep='\t')
write.table(umapClusters, file=paste0(home1Directory,'umapClusters_s1s2Scanorama_50pcs_filter_hg19.tsv'), col.names = TRUE, sep='\t')
write.table(logNormalizedCountsRound, file=paste0(home1Directory,'logNormalizedCounts_s1s2Scanorama_50pcs_filterRound_hg19.tsv'), col.names = TRUE, sep='\t')

rm(s1s2Scanorama)
rm(umapCoordinates)
rm(umapClusters)
rm(logNormalizedCounts)
rm(logNormalizedCountsRound)

#############################################################################
#DON'T NEED TO RUN BLOCK THIS AGAIN
#############################################################################