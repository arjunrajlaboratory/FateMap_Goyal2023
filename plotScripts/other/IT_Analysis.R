########################################################################################
#Description: Analysis of clusters from DEGs for Tirosh et al., 2016 data
#Author: Maalavika Pillai
#Version: 1
#Edited on: 7/11/22
########################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(svglite) 
library(ggrepel)
library(SeuratObject)
library(VennDiagram)
library(RColorBrewer)
library(tidyr)

dataDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/TiroshData/seurat/"
dataDirectory1 <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM01_DEGandVarGenes/Data/"
plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/TiroshData/seuratOutput/"


####Figure out : 
counts <- read.delim(paste0(dataDirectory,"GSE72056_melanoma_single_cell_revised_v2.txt"))
counts <- as.data.frame(t(counts))
names(counts) <- counts["Cell",]
cells <- rownames(counts)
counts <- apply(counts, 2, as.integer)
counts <- as.data.frame(counts)
rownames(counts) <- cells
#counts <- counts[counts$`malignant(1=no,2=yes,0=unresolved)`==2&counts$tumor==60,]
counts <- counts[counts$`malignant(1=no,2=yes,0=unresolved)`==2,]
counts <- counts[complete.cases(counts),]
tumor <- counts$tumor
counts <- counts[,colSums(counts)>10] #Remove genes with less than 10 counts
counts <- counts[,-c(1:2)]
counts <- as.data.frame(t(counts))

df <- CreateSeuratObject(counts = counts)
df <- NormalizeData(df, normalization.method = "LogNormalize")
df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 7000)
all.genes <- rownames(df)
df <- ScaleData(df)
df <- RunPCA(df)
df <- RunUMAP(df, dims = 1:50)
saveRDS(df, paste0(dataDirectory,"TiroshData_UMAP.rds"))

df <- readRDS(paste0(dataDirectory,"TiroshData_UMAP.rds"))
df <- FindNeighbors(df)
df <- FindClusters(df)
df.markers <- FindAllMarkers(df)
write.table(df.markers, paste0(dataDirectory, "TiroshData_clusterDEG.tsv") , sep = "\t", quote = F)
df.markers <- read.delim(paste0(dataDirectory, "TiroshData_clusterDEG.tsv") , sep = "\t",)
df.markers <-  df.markers %>% filter(p_val_adj <0.05& avg_log2FC>1)
df.markers<- unique(df.markers$gene)

clusterDEG_FM01 <- read.delim("~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM01_DEGandVarGenes/Data/scTransformMarkers_snn06.tsv", sep ="\t")
clusterDEG_FM01 <- clusterDEG_FM01 %>% filter(p_val_adj <0.05& avg_logFC>1)
clusterDEG_FM01 <- unique(clusterDEG_FM01$gene)
clusterDEG_FM01 <- clusterDEG_FM01[clusterDEG_FM01 %in% rownames(df@assays$RNA)]
  
#Fisher's exact test
fm01 <- readRDS("~/Library/CloudStorage/OneDrive-NorthwesternUniversity/Fatemap Revisions/Data/Drug treated and naive/FM01/s1s2ScTransform_50pcs_filter.rds")
nAllGenes = sum(rownames(df@assays$RNA) %in% rownames(fm01@assays$RNA))
nClusterDEG = length(clusterDEG_FM01)
nPatientDEG = length(df.markers)
nOverlap = length(intersect(clusterDEG_FM01,df.markers))
p <- sum(dhyper(nOverlap:nPatientDEG,nClusterDEG, nAllGenes-nClusterDEG, nPatientDEG)) #Hyper geometric dist: Probability that nOverlap or more genes are Overlapping in the two sets. 
p

#Venn Diagram
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(clusterDEG_FM01 ,df.markers),
  category.names = c( "Cell Line", "Patient" ),
  filename = paste0(plotDirectory, "ClusterDEG_CellLinevsPatient_Venn.png" ),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 550 , 
  width = 550 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  
  # Numbers
  cex = .6,
  fontface = "bold",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-23, 23),
  cat.dist = c(0.055, 0.055)
)

list = c("ACTG2", "ACTA2", "MYOCD", "TAGLN", "IFIT2", "OASL", "CCL3", "DDX58", "VCAM1", "PKDCC", "TDO2", "FOXF2", "NGFR",
         "COL9A3", "S100B", "ITGA6", "GAS7", "MLANA", "SOX10", "MITF", "PMEL", "AXL", "SERPINE1", "BGN", "WNT5A", "EGFR", "NRG1")
plot <- FeaturePlot(df,list[1:15])
ggsave(paste0(plotDirectory, "FeaturePlot1_TiroshData.svg"), plot, width = 20, height = 20)
plot <- FeaturePlot(df,list[-c(1:15)])
ggsave(paste0(plotDirectory, "FeaturePlot2_TiroshData.svg"), plot, width = 20, height = 15)


#Create random genes for each situation
df <- readRDS(paste0(dataDirectory,"TiroshData_UMAP.rds"))
df.markers <- read.delim(paste0(dataDirectory, "TiroshData_clusterDEG.tsv") , sep = "\t",)
df.markers <-  df.markers %>% filter(p_val_adj <0.05& avg_log2FC>1)
df.markers <- unique(df.markers$gene)
fm01 <- readRDS("~/Library/CloudStorage/OneDrive-NorthwesternUniversity/Fatemap Revisions/Data/Drug treated and naive/FM01/s1s2ScTransform_50pcs_filter.rds")
AllGenes <- intersect(rownames(df@assays$RNA@counts), rownames(fm01@assays$RNA@counts))
DEG_FM01 <- read.delim(paste0(dataDirectory1, "DEG_FC2_pval0.01.tsv"), sep = "\t")$gene
DEG_FM01 <- DEG_FM01[DEG_FM01 %in% AllGenes ]
fracOverlap_DEG <- length(intersect(df.markers, DEG_FM01))/ min(length(DEG_FM01), length(df.markers))
varGenes_FM01 <-  read.delim(paste0(dataDirectory1, "topVarGenes_1000.tsv"), sep = "\t")$x
varGenes_FM01 <- varGenes_FM01[varGenes_FM01 %in% AllGenes] 
fracOverlap_varGenes <- length(intersect(df.markers, varGenes_FM01))/ min(length(varGenes_FM01), length(df.markers))

fracOverlap <- data.frame()
for (i in 1:100) {
  set.seed(i)
  randGenes_overlap1 <- sample(AllGenes , 1000) 
  set.seed(100*i)
  randGenes_overlap2 <- sample(AllGenes , 1000) 
  fracOverlap[i,1:2] <- c(length(intersect(randGenes_overlap1, DEG_FM01))/min(length(randGenes_overlap1), length (varGenes_FM01)),
                          length(intersect(randGenes_overlap1, varGenes_FM01))/min(length(randGenes_overlap1), length (varGenes_FM01)))
}
names(fracOverlap) <- c("DEG", "VarGenes")
fracOverlap  <- fracOverlap %>% 
  pivot_longer(cols = DEG:VarGenes)
fracOverlap$sample <- "Random"
fracOverlap <- rbind(fracOverlap, data.frame("name" = "DEG", "value" = fracOverlap_DEG, "sample" = "Real"))
fracOverlap <- rbind(fracOverlap, data.frame("name" = "VarGenes", "value" = fracOverlap_varGenes , "sample" = "Real"))

set.seed(100)
plot <- ggplot(fracOverlap, aes(x = name, y = value, color = sample))+
  geom_jitter( position=position_jitterdodge(dodge.width =0.8, jitter.width = 0.2), shape = 16) +
  labs(x="", y = "Fraction overlap")+
  theme_classic() 
ggsave(plot = plot,filename= paste0(plotDirectory, "IT_Random_real_overlap.svg"), height=10, width=10)
write.csv(fracOverlap, paste0(plotDirectory, "IT_Random_real_overlap.csv"))



