########################################################################################
#Description: Cell cycle analysis for FM04
#Author: Maalavika Pillai
#Version: 1
#Edited on: 11/2/22
########################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(svglite) 
library(ggrepel)


dataDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM04/seurat/"
plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM04/seuratOutput/"

#Barcodes
######Barcode entries for Sample 1 and Sample 2
#################################################
barcode50 = as_tibble(read.table( '~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM04/barcodes/stepThreeStarcodeShavedReads.txt', stringsAsFactors=F, header = T))
names(barcode50)[names(barcode50) == "SampleNum"] <- "sampleNum" #Rename to keep it consistent with previously saved files
umiCut = 6 # Minimum UMI cutoff for reliable analysis. Can be played around with between 3 to 5.
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.  

upToLineageCounting = barcode50 %>% 
  filter(sampleNum %in% c("1","2")) %>%
  group_by(cellID, BC50StarcodeD8, sampleNum) %>%
  summarise(nUMI = length(sampleNum)) %>%
  filter(nUMI >= umiCut) %>% 
  group_by(cellID) %>%
  mutate(nLineages = length(cellID)) 

#making the nomenclature consistent with the Matrix from data integration
upToLineageCounting$sampleNum = sub("^", "S", upToLineageCounting$sampleNum) 
#####
###taking only single cellid-->barcode mappings
linCountTooverlaps = upToLineageCounting %>%
  ungroup() %>%
  filter(nLineages <= linCut) %>%
  unique()
barcodes  = upToLineageCounting %>%
  ungroup() %>%
  filter(nLineages <= linCut) %>%
  select(-nUMI,-nLineages) %>%
  unique()

write.table(barcodes, "barcodeCellID.tsv", quote = F, row.names = F, sep ="\t")
#########################33

#Load data and barcodes
fm04 <- readRDS(paste0(dataDirectory, "s1s2ScTransform_50pcs_filter.rds"))
barcodes <- read.delim(paste0(dataDirectory, "barcodeCellID.tsv"), sep = "\t")
barcodes$cellID <- paste0(barcodes$sampleNum, "_", barcodes$cellID)

#Extract cell cycle genes, from Tirosh et al., 2015
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Cell cycle scoring
fm04 <- CellCycleScoring(fm04, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#Check if UMAP localisation is dependent on cell cycle scores
plot <- DimPlot(fm04)
ggsave(paste0(plotDirectory, "FM04_Uncorrected_PhaseUMAP.svg"), width = 10, height = 10)

#Regress out cell cycle genes
fm04.cc <- ScaleData(fm04, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(fm04))
#Check if regression has occured by looking at PCA
fm04.cc <- RunPCA(fm04.cc, features = c(s.genes, g2m.genes))
DimPlot(fm04.cc, reduction = 'pca')
fm04.cc <- RunPCA(fm04.cc,  features = VariableFeatures(fm04))
fm04.cc <- RunUMAP(fm04.cc, dims = 1:50)
plot <- DimPlot(fm04.cc)
ggsave(paste0(plotDirectory, "FM04_CellCycleCorrected_PhaseUMAP.svg"), width = 10, height = 10)
saveRDS(fm04.cc, paste0(dataDirectory, "s1s2ScTransform_50pcs_filter_CellCycleCorrected.rds"))


#Look at dominant cluster fraction in each clone
####Domninant Cluster analysis

fm04.cc <- readRDS(file.choose())

fm04.cc <- FindNeighbors(fm04.cc)
fm04.cc <- FindClusters(object=fm04.cc, resolution = 0.5, verbose = FALSE)
umapCoordinates = (fm04.cc[['umap']])@cell.embeddings #UMAP coordinates
cells_UMAP = rownames(umapCoordinates) #CellIds with Sample number as prefix
cells_UMAP_cellID = sub("S\\d_", "", cells_UMAP)
cells_UMAP_sampleNum = gsub("[^S12]", "", cells_UMAP)
umapCoordinates = as_tibble(umapCoordinates) #UMAP coordinates in tibble
umapCoordinates = umapCoordinates %>% mutate(cellID = cells_UMAP_cellID,
                                             sampleNum = cells_UMAP_sampleNum)
jointUMAP <- inner_join(linCountTooverlaps, umapCoordinates, by = c("cellID", "sampleNum")) %>% select(-nLineages)

jointUMAPS1 = jointUMAP %>% filter(sampleNum == "S1")
jointUMAPS2 = jointUMAP %>% filter(sampleNum == "S2")

umapClusters = (fm04.cc[['seurat_clusters']])
cells_Clusters = rownames(umapClusters) #CellIds with Sample number as prefix
cells_Clusters_cellID = sub("S\\d_", "", cells_Clusters)
cells_Clusters_sampleNum = gsub("[^S12]", "", cells_Clusters)
umapClusters = as_tibble(umapClusters)
umapClusters = umapClusters %>% mutate(cellID = cells_Clusters_cellID,
                                       sampleNum = cells_Clusters_sampleNum)

umapClustersS1 = umapClusters %>% filter(sampleNum == "S2")
coloniesToAnalyze2 = jointUMAPS2 %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >4)
colonyUMAP = inner_join(coloniesToAnalyze2,jointUMAPS2, by = c("BC50StarcodeD8")) %>% select(BC50StarcodeD8,cellID,nColony)
colonyCluster = inner_join(colonyUMAP,umapClustersS1, by = c("cellID")) %>% select(-cellID,-sampleNum,nColony) %>% ungroup()

fractionColonies = colonyCluster %>% group_by(BC50StarcodeD8,seurat_clusters) %>% summarise(nCount = length(seurat_clusters))
fractionColonies1 <- fractionColonies %>% group_by(BC50StarcodeD8) %>% filter(nCount == max(nCount))%>% slice(1)
finalFractionColonies = inner_join(fractionColonies1,coloniesToAnalyze2, by = "BC50StarcodeD8") %>% mutate(nfractionFinal = nCount/nColony) %>% select(nfractionFinal) %>% mutate(type = "colonies")

fractionColoniesSeuratClusters = fractionColonies1 %>% select(BC50StarcodeD8,seurat_clusters)
fractionColoniesSeuratClustersColony = inner_join(fractionColoniesSeuratClusters,coloniesToAnalyze2, by = "BC50StarcodeD8")
#fractionColoniesA$seurat_clusters <- as.numeric(fractionColoniesA$seurat_clusters) # then check that fractionColoniesA looks the same as fractionColonies
#####
maxFraction = c()
medianFraction = c()
meanFraction = c()


for (i in c(1:length(fractionColoniesSeuratClustersColony$BC50StarcodeD8))) {
  for (j in c(1:50)) {
    subSampled = sample_n(colonyCluster, fractionColoniesSeuratClustersColony$nColony[i])
    fraction = subSampled %>% group_by(seurat_clusters) %>% summarise(nFraction = length(seurat_clusters)/fractionColoniesSeuratClustersColony$nColony[i])
    if (fractionColoniesSeuratClustersColony$seurat_clusters[i] %in% fraction$seurat_clusters) {
      fractionTemp = fraction %>% filter(seurat_clusters == fractionColoniesSeuratClustersColony$seurat_clusters[i])
      maxFraction[j] = fractionTemp$nFraction
    } else {
      maxFraction[j] = 0
    }
  }
  meanFraction[i] = mean(maxFraction)
  medianFraction[i] = median(maxFraction) 
}

finalFractionRandom = tibble(nfractionFinal = numeric()) %>% add_row(nfractionFinal = meanFraction) %>% mutate(type = "random")
fractionFinalClusters = bind_rows(finalFractionColonies,finalFractionRandom)
fractionFinalClusters$type <- factor(fractionFinalClusters$type, levels=c("random", "colonies"))
write.table(fractionFinalClusters, file=paste0(dataDirectory,'FM04_fractionFinalClusters_S2_snn05_CellCycleCorrected.tsv'), col.names = TRUE, sep='\t')

summaryfractionFinalClusters = fractionFinalClusters %>% group_by(type) %>% summarise(mean = mean(nfractionFinal))

#### for S1, random = 0.0919; observed = 0.387
plot = ggplot(fractionFinalClusters, aes(x = type, y=nfractionFinal, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  geom_jitter(width = 0.1, shape = 16) +
  stat_compare_means( size =4, label.x.npc = "left")+
  scale_fill_manual(values=c("grey","goldenrod3")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none") +
  ylim(0,1)
fractionFinalClusters %>% group_by(type) %>% summarize(mean(nfractionFinal))
ggsave(plot, file = paste0(plotDirectory, 'FM04_ColonySpecificPlotColony_V1_withDotsFinal_S2_05_ CellCycleCorrected.svg'), width = 4, height = 8)



#####3 largest clones
largestBarcodes <- coloniesToAnalyze2$BC50StarcodeD8[order(coloniesToAnalyze2$nColony, decreasing = T)]
largestBarcodes <- largestBarcodes[1:3]
plots <- list()
for (i in 1:length(largestBarcodes)) {
  cell_names <- paste0("S2_",colonyUMAP$cellID[colonyUMAP$BC50StarcodeD8 %in% largestBarcodes[i]])
  plots[[i]] <- DimPlot(fm04.cc, cells.highlight = cell_names , cols.highlight = "hotpink3") + NoLegend()
  
}
plot <- ggarrange(plotlist = plots, nrow = 1)
ggsave(plot = plot, filename = paste0(plotDirectory,"FM04_UMAP_specificClones.svg"), width = 12, height = 3)
