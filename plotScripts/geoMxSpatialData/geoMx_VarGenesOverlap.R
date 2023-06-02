########################################################################################
#Description: Check overlap of variable genes in spatial transcriptomics data 
#and FM01
#Author: Maalavika Pillai
#Version: 1
#Edited on: 8th December 2022
########################################################################################
library(dplyr)
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

VennDiag <- function(geneList1, geneList2,name1, name2, plotDirectory){
  myCol <- brewer.pal(3, "Pastel2")
  venn.diagram(
    x = list(geneList1 ,geneList2),
    category.names = c( name1,name2 ),
    filename = paste0(plotDirectory, "Venn_", name1,"_",name2, ".png" ),
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
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-23, 23),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans"
  )
  
}

dataDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/geoMX Data/"
dataDirectory1 <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM01_DEGandVarGenes/Data/"
plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/geoMX Data/Plots/"

varGenes_geoMx <- read.csv(paste0(dataDirectory, "geomx_variable_genes_by_patient_condition.csv"))

#Read pre-existing DEGs and variable genes
DEG_FM01 <- read.delim(paste0(dataDirectory1, "DEG_FC2_pval0.01.tsv"), sep = "\t")$gene
DEG_FM01 <- DEG_FM01[DEG_FM01 %in% varGenes_geoMx$gene]
varGenes_FM01 <-  read.delim(paste0(dataDirectory1, "topVarGenes_1000.tsv"), sep = "\t")$x
varGenes_FM01 <- varGenes_FM01[varGenes_FM01 %in% varGenes_geoMx$gene] 

#For fisher's exact test
fm01 <- readRDS("~/Library/CloudStorage/OneDrive-NorthwesternUniversity/Fatemap Revisions/Data/Drug treated and naive/FM01/s1s2ScTransform_50pcs_filter.rds")
nAllGenes_GeoMx <- length(unique(varGenes_geoMx$gene))
nAllGenes = length(intersect(unique(varGenes_geoMx$gene),rownames(fm01@assays$RNA)))
nClusterDEG = length(DEG_FM01)
nVarGenes = length(varGenes_FM01)
nGeoMxVarGenes = 1000

#Get list of DEGs / variable genes that are present in the geoMx data
varGenes_geoMx <- varGenes_geoMx %>% filter(gene %in% rownames(fm01@assays$RNA))
varGenes_geoMx$surg_immune_pre <- gsub("SURG_IMMUNE", "SURG",varGenes_geoMx$surg_immune_pre)
geoMx_og <- varGenes_geoMx 

varGenes_geoMx <- varGenes_geoMx %>%
  group_by(patient_id, surg_immune_pre) %>%
  slice_max(cv_cpm, n=1000)


varGenes_overlap <- varGenes_geoMx %>% 
  group_by(patient_id, surg_immune_pre) %>%
  summarize("overlapFraction" = length(intersect(gene, varGenes_FM01))/ min(length(gene), length (varGenes_FM01)), "overlapCount" =  length(intersect(gene, varGenes_FM01)))

DEGGenes_overlap <- varGenes_geoMx %>% 
  group_by(patient_id, surg_immune_pre) %>%
  summarize("overlapFraction" = length(intersect(gene, DEG_FM01))/  min(length(gene), length (DEG_FM01)), "overlapCount" =  length(intersect(gene, DEG_FM01)))

#Create Venn diagrams
#For surg wrt FM01
varGenes_geoMx_surg <- varGenes_geoMx %>%
  filter(surg_immune_pre== "SURG")
p_DEG <- vector()
p_Var <- vector()
for (i in unique(varGenes_geoMx_surg$patient_id)) {
    genes <- varGenes_geoMx_surg$gene[varGenes_geoMx_surg$patient_id==i]
    VennDiag(genes, DEG_FM01, "geoMx", "DEG_CellLine", paste0(plotDirectory,"geoMx_", i, "_"))
    nOverlap = length(intersect(genes, DEG_FM01))
    p_DEG[paste0(i)] <- sum(dhyper(nOverlap:nGeoMxVarGenes,nClusterDEG, nAllGenes-nClusterDEG, nGeoMxVarGenes)) #Hyper geometric dist: Probability that nOverlap or more genes are Overlapping in the two sets. 
    
    VennDiag(genes, varGenes_FM01, "geoMx", "VarGenes_CellLine", paste0(plotDirectory,"geoMx_", i, "_"))
    nOverlap = length(intersect(genes, varGenes_FM01))
    p_Var[paste0(i)]  <- sum(dhyper(nOverlap:nGeoMxVarGenes, nVarGenes, nAllGenes-nVarGenes, nGeoMxVarGenes)) #Hyper geometric dist: Probability that nOverlap or more genes are Overlapping in the two sets. 
}
write.csv(cbind(p_DEG,p_Var), paste0(plotDirectory,"geoMx_p_Value_Hypergeometric_FM01_vs_geoMx.csv"))

#For pre wrt sug
varGenes_geoMx_pre_surg <- varGenes_geoMx %>%
  filter(patient_id %in% c(109,134))
p_pre_surg <- vector()
varGenes_geoMx_pre_surg_overlap <- varGenes_geoMx_pre_surg %>%
  group_by(patient_id) %>%
  summarize("overlapFraction" = length(intersect(gene[surg_immune_pre == "PRE"],gene[surg_immune_pre == "SURG"]))/ 1000,
            "overlapCount" =  length(intersect(gene[surg_immune_pre == "PRE"],gene[surg_immune_pre == "SURG"])))
for (i in unique(varGenes_geoMx_pre_surg$patient_id)) {
  VennDiag(varGenes_geoMx_pre_surg$gene[varGenes_geoMx_pre_surg$patient_id==i&varGenes_geoMx_pre_surg$surg_immune_pre=="PRE"], 
           varGenes_geoMx_pre_surg$gene[varGenes_geoMx_pre_surg$patient_id==i&varGenes_geoMx_pre_surg$surg_immune_pre=="SURG"],
           "geoMx_pre", "geoMx_surg", paste0(plotDirectory,"geoMx_", i, "_"))
  nOverlap <- length(intersect(varGenes_geoMx_pre_surg$gene[varGenes_geoMx_pre_surg$patient_id==i&varGenes_geoMx_pre_surg$surg_immune_pre=="PRE"], 
                        varGenes_geoMx_pre_surg$gene[varGenes_geoMx_pre_surg$patient_id==i&varGenes_geoMx_pre_surg$surg_immune_pre=="SURG"]))
  
  p_pre_surg[as.character(i)] <- sum(dhyper(nOverlap:1000, 1000, (nAllGenes_GeoMx-1000), 1000)) 
}
write.csv(p_pre_surg, paste0(plotDirectory,"geoMx_p_Value_Hypergeometric_prevspost.csv"))

#Create random genes for each situation
fracOverlap <- data.frame()
for (i in 1:100) {
  set.seed(i)
  randGenes_overlap1 <- sample(geoMx_og$gene, 1000) 
  set.seed(100*i)
  randGenes_overlap2 <- sample(geoMx_og$gene, 1000) 
  fracOverlap[i,1:3] <- c(length(intersect(randGenes_overlap1, DEG_FM01))/min(length(randGenes_overlap1), length (varGenes_FM01)),length(intersect(randGenes_overlap1, varGenes_FM01))/min(length(randGenes_overlap1), length (varGenes_FM01)), 
                         length(intersect(randGenes_overlap1, randGenes_overlap2))/1000)
}
names(fracOverlap) <- c("DEG", "VarGenes", "Pre_Post")
fracOverlap  <- fracOverlap %>% 
  pivot_longer(cols = DEG:Pre_Post)
fracOverlap$sample <- "Random"
fracOverlap <- rbind(fracOverlap, data.frame("name" = "DEG", "value" = DEGGenes_overlap$overlapFraction[DEGGenes_overlap$surg_immune_pre=="SURG"], "sample" = "Real"))
fracOverlap <- rbind(fracOverlap, data.frame("name" = "VarGenes", "value" = varGenes_overlap$overlapFraction[varGenes_overlap$surg_immune_pre=="SURG"], "sample" = "Real"))
fracOverlap <- rbind(fracOverlap, data.frame("name" = "Pre_Post", "value" = varGenes_geoMx_pre_surg_overlap$overlapFraction, "sample" = "Real"))

fracOverlap_pre_post <- fracOverlap %>%
  filter(name =="Pre_Post" )
mean_line <-fracOverlap_pre_post %>% group_by(sample) %>% summarise(mean_y = mean(value));
set.seed(100)
plot <- ggplot(fracOverlap_pre_post, aes(x = sample, y = value, color = sample))+
  geom_jitter( position=position_jitterdodge(dodge.width =0.8, jitter.width = 0.2), shape = 16) +
  labs(x="", y = "Fraction overlap")+
  geom_path(data = mean_line, aes(x = sample, y = mean_y, group = 1), color = "grey")+
  theme_classic() 
ggsave(plot = plot,filename= paste0(plotDirectory, "geoMx_Random_real_prePost_overlap.svg"), height=20, width=8)

fracOverlap <- fracOverlap %>%
  filter(name =="VarGenes" )
mean_line <-fracOverlap %>% group_by(sample) %>% summarise(mean_y = mean(value));
set.seed(100)
plot <- ggplot(fracOverlap, aes(x = sample, y = value, color = sample))+
  geom_jitter( position=position_jitterdodge(dodge.width =0.8, jitter.width = 0.2), shape = 16) +
  labs(x="", y = "Fraction overlap")+
  geom_path(data = mean_line, aes(x = sample, y = mean_y, group = 1), color = "grey")+
  theme_classic() 
ggsave(plot = plot,filename= paste0(plotDirectory, "geoMx_Random_real_overlap.svg"), height=20, width=8)
  
  