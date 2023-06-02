####E-distance

library(scperturbR)
library(Seurat)
subsamplingNumber <- 10
resultPath <- '~/Dropbox (RajLab)/FateMap/FINAL_REVISION/analysis_Jonas/data/e_dist'

#load data
sampleFM011.data <- Read10X(data.dir = "~/Dropbox (RajLab)/FateMap/FINAL_REVISION/analysis_Jonas/rawData/FM01/1_1uMPLX/filtered_feature_bc_matrix/")
sampleFM012.data <- Read10X(data.dir = "~/Dropbox (RajLab)/FateMap/FINAL_REVISION/analysis_Jonas/rawData/FM01/2_1uMPLX/filtered_feature_bc_matrix/")
sampleFM061.data <- Read10X(data.dir = "~/Dropbox (RajLab)/FateMap/FINAL_REVISION/analysis_Jonas/rawData/FM06/WM989Naive_1/filtered_feature_bc_matrix/")
sampleFM062.data <- Read10X(data.dir = "~/Dropbox (RajLab)/FateMap/FINAL_REVISION/analysis_Jonas/rawData/FM06/WM989Naive_2/filtered_feature_bc_matrix/")

#Preprocess and create Seurat objects
fm011.seurat <- CreateSeuratObject(counts = sampleFM011.data, project = "FM01_1_1uMPLX", min.cells = 3, min.features = 200)
fm011.seurat[["percent.mt"]] <- PercentageFeatureSet(object = fm011.seurat, pattern = "^MT-")
fm011.seurat <- subset(x = fm011.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)

fm012.seurat <- CreateSeuratObject(counts = sampleFM012.data, project = "FM01_2_1uMPLX", min.cells = 3, min.features = 200)
fm012.seurat[["percent.mt"]] <- PercentageFeatureSet(object = fm012.seurat, pattern = "^MT-")
fm012.seurat <- subset(x = fm012.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 26)

fm061.seurat <- CreateSeuratObject(counts = sampleFM061.data, project = "WM989Naive_1", min.cells = 3, min.features = 200)
fm061.seurat[["percent.mt"]] <- PercentageFeatureSet(object = fm061.seurat, pattern = "^MT-")
fm061.seurat <- subset(x = fm061.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 20)

fm062.seurat <- CreateSeuratObject(counts = sampleFM062.data, project = "WM989Naive_2", min.cells = 3, min.features = 200)
fm062.seurat[["percent.mt"]] <- PercentageFeatureSet(object = fm062.seurat, pattern = "^MT-")
fm062.seurat <- subset(x = fm062.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7200 & percent.mt < 20)

#get minimum and maximum number of cells per sample
ncells <- min(ncol(fm011.seurat), ncol(fm012.seurat),ncol(fm061.seurat), ncol(fm062.seurat) )
maxCells <- max(ncol(fm011.seurat), ncol(fm012.seurat),ncol(fm061.seurat), ncol(fm062.seurat))

#calculate and save e-distances for different resolutions
analysis.seurat <- fm061.seurat
set.seed(123)
random <- sample.int(1000, subsamplingNumber)
for (s in 1:subsamplingNumber)
{
  set.seed(random[s])
  cells <- sample(colnames(analysis.seurat), size =ncells, replace=F)
  analysis.seurat.subset <- analysis.seurat[, cells]
  analysis.seurat.subset <- SCTransform(analysis.seurat.subset, verbose = FALSE)
  analysis.seurat.subset <- RunPCA(object = analysis.seurat.subset, assay = "SCT")
  analysis.seurat.subset <- FindNeighbors(object=analysis.seurat.subset, dims=1:50, reduction = "pca", force.recalc = TRUE)
  for (i in 1:10)
  {
    resolution = 0.1*i
    analysis.seurat.subset <- FindClusters(object=analysis.seurat.subset, resolution=resolution, verbose = FALSE)
    analysis.seurat.subset <- RunUMAP(object = analysis.seurat.subset, reduction = "pca", dims = 1:50, reduction.name = "umap", verbose = FALSE)
    estats <- edist(analysis.seurat.subset, groupby = paste0("SCT_snn_res.", as.character(resolution)))
    dir.create(file.path(resultPath, "fm06_1"), showWarnings = FALSE)
    write.csv(estats, paste0(resultPath, "/fm06_1/","Edist_fm06_1_res", as.character(resolution), "_", as.character(s)  ,".csv"))
  }
}

analysis.seurat <- fm011.seurat
set.seed(123)
random <- sample.int(1000, subsamplingNumber)
for (s in 1:subsamplingNumber)
{
  set.seed(random[s])
  cells <- sample(colnames(analysis.seurat), size =ncells, replace=F)
  analysis.seurat.subset <- analysis.seurat[, cells]
  analysis.seurat.subset <- SCTransform(analysis.seurat.subset, verbose = FALSE)
  analysis.seurat.subset <- RunPCA(object = analysis.seurat.subset, assay = "SCT")
  analysis.seurat.subset <- FindNeighbors(object=analysis.seurat.subset, dims=1:50, reduction = "pca", force.recalc = TRUE)
  for (i in 1:10)
  {
    resolution = 0.1*i
    analysis.seurat.subset <- FindClusters(object=analysis.seurat.subset, resolution=resolution, verbose = FALSE)
    analysis.seurat.subset <- RunUMAP(object = analysis.seurat.subset, reduction = "pca", dims = 1:50, reduction.name = "umap", verbose = FALSE)
    estats <- edist(analysis.seurat.subset, groupby = paste0("SCT_snn_res.", as.character(resolution)))
    dir.create(file.path(resultPath, "fm01_1"), showWarnings = FALSE)
    write.csv(estats, paste0(resultPath, "/fm01_1/","Edist_fm01_1_res", as.character(resolution), "_", as.character(s)  ,".csv"))
  }
}

analysis.seurat <- fm012.seurat
set.seed(123)
random <- sample.int(1000, subsamplingNumber)
for (s in 1:subsamplingNumber)
{
  set.seed(random[s])
  cells <- sample(colnames(analysis.seurat), size =ncells, replace=F)
  analysis.seurat.subset <- analysis.seurat[, cells]
  analysis.seurat.subset <- SCTransform(analysis.seurat.subset, verbose = FALSE)
  analysis.seurat.subset <- RunPCA(object = analysis.seurat.subset, assay = "SCT")
  analysis.seurat.subset <- FindNeighbors(object=analysis.seurat.subset, dims=1:50, reduction = "pca", force.recalc = TRUE)
  for (i in 1:10)
  {
    resolution = 0.1*i
    analysis.seurat.subset <- FindClusters(object=analysis.seurat.subset, resolution=resolution, verbose = FALSE)
    analysis.seurat.subset <- RunUMAP(object = analysis.seurat.subset, reduction = "pca", dims = 1:50, reduction.name = "umap", verbose = FALSE)
    estats <- edist(analysis.seurat.subset, groupby = paste0("SCT_snn_res.", as.character(resolution)))
    dir.create(file.path(resultPath, "fm01_2"), showWarnings = FALSE)
    write.csv(estats, paste0(resultPath, "/fm01_2/","Edist_fm01_2_res", as.character(resolution), "_", as.character(s)  ,".csv"))
  }
}

analysis.seurat <- fm062.seurat
set.seed(123)
random <- sample.int(1000, subsamplingNumber)
for (s in 1:subsamplingNumber)
{
  set.seed(random[s])
  cells <- sample(colnames(analysis.seurat), size =ncells, replace=F)
  analysis.seurat.subset <- analysis.seurat[, cells]
  analysis.seurat.subset <- SCTransform(analysis.seurat.subset, verbose = FALSE)
  analysis.seurat.subset <- RunPCA(object = analysis.seurat.subset, assay = "SCT")
  analysis.seurat.subset <- FindNeighbors(object=analysis.seurat.subset, dims=1:50, reduction = "pca", force.recalc = TRUE)
  for (i in 1:10)
  {
    resolution = 0.1*i
    analysis.seurat.subset <- FindClusters(object=analysis.seurat.subset, resolution=resolution, verbose = FALSE)
    analysis.seurat.subset <- RunUMAP(object = analysis.seurat.subset, reduction = "pca", dims = 1:50, reduction.name = "umap", verbose = FALSE)
    estats <- edist(analysis.seurat.subset, groupby = paste0("SCT_snn_res.", as.character(resolution)))
    dir.create(file.path(resultPath, "fm06_2"), showWarnings = FALSE)
    write.csv(estats, paste0(resultPath, "/fm06_2/","Edist_fm06_2_res", as.character(resolution), "_", as.character(s)  ,".csv"))
  }
}

.
