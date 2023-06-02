########Adapted from Correlation_clone_vs_random_all.R (orginally for FM01) 
#and FM08_varGenes.R for FM10

########################################################################################
#Description: Homogeneity analysis for FateMap data
#Aim: Estimate pairwise correlation and compare for clones vs non-clones
#Author: Maalavika Pillai
#Version: 5.2
#Created on: 11/2/22
#Use top n variable genes and all genes, across different clone sizes
########################################################################################

library(Seurat)
library(ggpubr)
library(dplyr)
library(rstatix)
library(ggplot2)
library(tidyr)

#Load data and barcodes
dataDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM10/jointAnalysis/seurat/"
plotDirectory <- "~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2022_FM10/jointAnalysis/seuratOutput/"

#Load data and barcodes
fm10 <- readRDS(paste0(dataDirectory, "s1s2ScTransform_50pcs_filter.rds"))
barcodes <- read.table(paste0(dataDirectory,'barcodeCellID.tsv'), stringsAsFactors=F, header = T, sep="\t")
barcodes$cellID <- paste0(barcodes$sampleNum, "_", barcodes$cellID)
fm10 <- RenameCells(fm10, new.names = substr(colnames(fm10),1,nchar(colnames(fm10))-2))

#Select samples that have barcode information available
df <- as.data.frame(t(as.matrix(fm10@assays$integrated@scale.data)))
 

#Remove barcodes with single cellIDs and only those absent in df
barcodes <- barcodes %>%
  filter(cellID %in% rownames(df)) %>%
  group_by(BC50StarcodeD8) %>%
  filter(length(cellID) >= 2)
df <- df[barcodes$cellID,]
barcodes$BC50StarcodeD8 <- factor(barcodes$BC50StarcodeD8, labels = 1:length(unique(barcodes$BC50StarcodeD8)))

#Calculate correlation matrix for all samples (useful if you plan on running iterations of subsampling)
corr <- cor(t(df), method="spearman")

#Calculate mean corrleation across pairs of cells within a clone 
#Paired control is selected as equal number of cells
#Use different clone size cut-offs
clone_corr_mat <- data.frame()

for (clone_size in c(2,5,10,20,50,100)) { #Max size <20
#for (clone_size in c(2,5,10)) {
  corr_mat <- data.frame()
  nBarcodes <- levels(barcodes$BC50StarcodeD8)[which(table(barcodes$BC50StarcodeD8)>clone_size)]
  set.seed(10)
  seed <- sample(1:10000, length(nBarcodes)) #>1000
  for(i in nBarcodes){
    pairs <- barcodes$cellID[barcodes$BC50StarcodeD8==i]  #Select all cells in clone
    set.seed(seed[which(nBarcodes==i)])
    random <- sample(barcodes$cellID, sum(barcodes$BC50StarcodeD8==i)) #Select same number of cells randomly
    corr_pairs <-  mean(corr[pairs,pairs]) #Get correlation of all pairwise combinations within clone
    corr_random <- mean(corr[random,random]) #Get correlation of all pairwise combinations across random samples
    k=which(unique(barcodes$BC50StarcodeD8)==i)
    corr_mat[((2*k)-1),1:3] <- c(mean(corr_pairs), "Clone",k)
    corr_mat[(2*k),1:3] <- c(mean(corr_random), "Random",k)
  }
  clone_corr_mat <- rbind(clone_corr_mat, data.frame(corr_mat, "CloneSize"=clone_size))
}

clone_corr_mat <- clone_corr_mat[complete.cases(clone_corr_mat),]
names(clone_corr_mat)[1:3] <- c("Correlation", "Clonal_information", "Index")
clone_corr_mat$Clonal_information <- as.factor(clone_corr_mat$Clonal_information)
clone_corr_mat$Correlation <- as.numeric(clone_corr_mat$Correlation)


#Plot the data
#All clone sizes

plot <- ggplot(clone_corr_mat, aes( x = Clonal_information, y = Correlation, color = Clonal_information)) + 
  geom_jitter(alpha=0.8, width = 0.2)+
  stat_compare_means(paired = TRUE,tip.length = 0.02, size = 3)+
  scale_color_manual(values = c(blues9[4],blues9[8]))+
  theme_classic() + 
  stat_summary(fun.y= mean, geom="crossbar", width=0.5, color="black")+
  geom_line(aes(group=Index), alpha =0.3, color = "darkgray")+
  stat_summary(fun.y=mean, colour="black", alpha =0.6, geom="line", aes(group = 1))+
  facet_wrap(~CloneSize)+
  rremove("legend") + 
  xlab("")
ggsave(paste0(plotDirectory,"FM10_CorrelationClonesvsRandom_allGenes.svg"),plot)


#Foldchange
fc <- group_by(clone_corr_mat, CloneSize, Index) %>%
  summarise(FoldChange = Correlation[Clonal_information=="Clone"]/Correlation[Clonal_information=="Random"])

fc <- group_by(fc, CloneSize) %>%
  summarize(mean(FoldChange))
write.table(fc,paste0(plotDirectory,"FM10_CorrelationFC_allGenes.tsv"))

##########################################################

##########################################################
#For top n variable genes
#specify genes
df_main = df
for (nVarGenes in c(50,100,500, 1000)) {
  #specify genes
  fm10 <- FindVariableFeatures(fm10, assay="RNA", nfeatures = nVarGenes)
  top10 <- head(VariableFeatures(fm10, assay="RNA"), 30)
  plot1 <- VariableFeaturePlot(fm10, assay="RNA")
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  ggsave(paste0(plotDirectory,"FM10_VariableGenes_",nVarGenes,".svg"),plot2, width = 8, height= 10)
  
  VarGenes = VariableFeatures(fm10, assay="RNA")
  df= df_main
  df <- df[,names(df) %in% VarGenes]
  #Calculate correlation matrix for all samples (useful if you plan on running iterations of subsampling)
  corr <- cor(t(df), method="spearman")
  
  #Calculate mean corrleation across pairs of cells within a clone 
  #Paired control is selected as equal number of cells
  #Use different clone size cut-offs
  clone_corr_mat <- data.frame()
  
  for (clone_size in c(2,5,10)) {
    corr_mat <- data.frame()
    nBarcodes <- levels(barcodes$BC50StarcodeD8)[which(table(barcodes$BC50StarcodeD8)>clone_size)] 
    set.seed(10)
    seed <- sample(1:10000, length(nBarcodes))
    for(i in nBarcodes){
      pairs <- barcodes$cellID[barcodes$BC50StarcodeD8==i]  #Select all cells in clone
      set.seed(seed[which(nBarcodes == i)])
      random <- sample(barcodes$cellID, sum(barcodes$BC50StarcodeD8==i)) #Select same number of cells randomly
      corr_pairs <-  mean(corr[pairs,pairs]) #Get correlation of all pairwise combinations within clone
      corr_random <- mean(corr[random,random]) #Get correlation of all pairwise combinations across random samples
      k=which(unique(barcodes$BC50StarcodeD8)==i)
      corr_mat[((2*k)-1),1:3] <- c(mean(corr_pairs), "Clone",k)
      corr_mat[(2*k),1:3] <- c(mean(corr_random), "Random",k)
    }
    clone_corr_mat <- rbind(clone_corr_mat, data.frame(corr_mat, "CloneSize"=clone_size))
  }
  
  clone_corr_mat <- clone_corr_mat[complete.cases(clone_corr_mat),]
  names(clone_corr_mat)[1:3] <- c("Correlation", "Clonal_information", "Index")
  clone_corr_mat$Clonal_information <- as.factor(clone_corr_mat$Clonal_information)
  clone_corr_mat$Correlation <- as.numeric(clone_corr_mat$Correlation)
  
  
  #Plot the data
  #All clone sizes
  
  plot <- ggplot(clone_corr_mat, aes( x = Clonal_information, y = Correlation, color = Clonal_information)) + 
    geom_jitter(alpha=0.8, width = 0.2)+
    #stat_pvalue_manual(stat.test,tip.length = 0.02,bracket.nudge.y = 0.1, size = 8)+
    stat_compare_means(paired = TRUE,tip.length = 0.02, size = 3)+
    scale_color_manual(values = c(blues9[4],blues9[8]))+
    theme_classic() + 
    stat_summary(fun.y= mean, geom="crossbar", width=0.5, color="black")+
    #geom_boxplot(alpha=0, outlier.shape = NA, position = position_dodge(1), width = 0.8)+
    #geom_line(aes(group=Index), alpha =0.3, color = "darkgray")+
    stat_summary(fun.y=mean, colour="black", alpha =0.6, geom="line", aes(group = 1))+
    facet_wrap(~CloneSize)+
    rremove("legend") + 
    xlab("")
  ggsave(paste0(plotDirectory,"FM10_CorrelationClonesvsRandom_VarGenes",nVarGenes,".svg"),plot, width =6, height = 5)
  
  
  #Foldchange
  fc <- group_by(clone_corr_mat, CloneSize, Index) %>%
    summarise(FoldChange = Correlation[Clonal_information=="Clone"]/Correlation[Clonal_information=="Random"])
  
  fc <- group_by(fc, CloneSize) %>%
    summarize(mean(FoldChange))
  write.table(x = fc, file = paste0(plotDirectory,"FM10_CorrelationFC_VarGenes",nVarGenes,".tsv"),quote = F)
}

##########################################################
