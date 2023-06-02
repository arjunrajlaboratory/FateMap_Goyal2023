########################################################################################
#Description: SVM and differentially expressed genes
#Aim: Identify DEGs (2 v/s n clusters) and train SVM to extract clonal information
#Author: Maalavika Pillai
#Version: 3
#Created on: 11th January 2023
########################################################################################

library(Seurat)
library(ggpubr)
library(dplyr)
library(e1071)
library(ggcorrplot)
library(svglite)

#function to extract training barcodes
trainingSubset <- function(barcodes, i = 1){
  x <- sample(barcodes$cellID, floor(0.8*nrow(barcodes)))
  if(all(unique(barcodes$BC50StarcodeD8) %in% barcodes$BC50StarcodeD8[barcodes$cellID %in% x])){
    return(x)
  }else{
    i = i+1
    print(i)
    trainingSubset(barcodes)
  }
}

#Load SeuratObject adn barcodes
setwd("~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM01/Data")
fm01 <- readRDS("s1s2ScTransform_50pcs_filter.rds")
barcodes <- read.table('barcodeCellID.tsv', stringsAsFactors=F, header = T, sep="\t")
#For S2
#barcodes <- read.table('barcodeCellID_S2.tsv', stringsAsFactors=F, header = T, sep="\t")

barcodes$cellID <- paste0(barcodes$sampleNum, "_", barcodes$cellID)
barcodes <- barcodes %>%
  group_by(BC50StarcodeD8) %>%
  mutate(nSize = length(cellID))
barcodes <- barcodes[barcodes$cellID %in% colnames(fm01),]
fm01$barcodes <- barcodes$BC50StarcodeD8[match(rownames(fm01@meta.data),barcodes$cellID)]
#barcodes <- barcodes[barcodes$nSize >10,]
df <- as.data.frame(t(as.matrix(fm01@assays$integrated@scale.data)))
df$barcodes <- as.factor(fm01$barcodes)
df <- df[complete.cases(df),]

#Read DEGs 
datadir <- "/Users/meep/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Maalavika/FM01_DEGandVarGenes/Data/"
DEG <- read.delim(paste0(datadir, "CloneGenes_FC1_p0.05.tsv"), sep = "\t") #log2FC = 1, p < 0.05
DEG <- unique(DEG$gene)


#Top 10 clones
barcodes_sub <- barcodes[barcodes$nSize>100,]
set.seed(100)
train <- trainingSubset(barcodes_sub)
test <- barcodes_sub$cellID[!barcodes_sub$cellID %in% train]
fm01_train <- df[train,c(DEG, "barcodes")]
fm01_test <- df[test,DEG]
ground_truth <- df[test, "barcodes"]

#Train SVM
model <- svm(barcodes~. , data = fm01_train, type = "C")
prediction <- predict(model, fm01_test)
prediction <- factor(prediction, levels = unique(ground_truth), labels = 1:length(unique(ground_truth)))
ground_truth <- factor(ground_truth, levels = unique(ground_truth), labels = 1:length(unique(ground_truth)))
confusionMat <- table(as.character(prediction), as.character(ground_truth))
confusionMat <- apply(confusionMat, 2, function(x){100*x/sum(x)})
plot <- ggcorrplot(confusionMat, lab = T) + 
        labs(y= "Expected clone", x= "Observed clone")+ 
        theme(
          axis.title.x = element_text(angle = 0, color = 'grey20', size = 22),
          axis.title.y = element_text(angle = 90, color = 'grey20', size = 22)
        )+
        scale_fill_gradient2(limit = c(0,100), low = "white", high =  blues9[6], mid = blues9[3], midpoint = 50)
ggsave(plot = plot, filename = "../../../Plots/FM01/SVMconfusionMatrix_clones.svg", width = 5,height = 5)
#For S2
#ggsave(plot = plot, filename = "../Plots/SVMconfusionMatrix_clones_S2.svg", width = 5,height = 5)

#Compare to random
fm01_train_random <- df[train,c(DEG, "barcodes")]
fm01_train_random$barcodes <- sample(fm01_train_random$barcodes)
ground_truth <- df[test, "barcodes"]

#Train SVM
set.seed(50)
model_random <- svm(barcodes~. , data = fm01_train_random, type = "C")
prediction_random <- predict(model_random, fm01_test)
prediction_random<- factor(prediction_random, levels = unique(ground_truth), labels = 1:5)
ground_truth <- factor(ground_truth, levels = unique(ground_truth), labels = 1:5)
confusionMat_random <- table(as.character(prediction_random), as.character(ground_truth))
confusionMat_random <- apply(confusionMat_random, 2, function(x){100*x/sum(x)})
missing <- setdiff(c(1:ncol(confusionMat_random)),rownames(confusionMat_random))
if(length(missing) > 0){
  missing <- data.frame(matrix(rep(0,(5*length(missing))), nrow = length(missing)), row.names = missing)
  names(missing) <- c(1:ncol(confusionMat_random))
  confusionMat_random <- rbind(confusionMat_random,missing )
}
plot <- ggcorrplot(confusionMat_random, lab = T) + 
        labs(y= "Expected clone", x= "Observed clone")+ 
        theme(
          axis.title.x = element_text(angle = 0, color = 'grey20', size = 22),
          axis.title.y = element_text(angle = 90, color = 'grey20', size = 22)
        )+
        scale_fill_gradient2(limit = c(0,100), low = "white", high =  blues9[6], mid = blues9[3], midpoint = 50)
ggsave(plot = plot, filename = "../Plots/SVMconfusionMatrix_random.svg", width = 5,height = 5)
#For S2
#ggsave(plot = plot, filename = "../Plots/SVMconfusionMatrix_random_S2.svg", width = 5,height = 5)

#Jitter plot comparing the two
accuracy <- data.frame("Accuracy" = c(diag(as.matrix(confusionMat)), diag(as.matrix(confusionMat_random))), "Sample" = c(rep("Clones", ncol(confusionMat)), rep("Shuffled", ncol(confusionMat_random))))
plot <- ggplot(accuracy, mapping = aes(x = Sample, y = Accuracy, color = Sample))+
  geom_jitter(alpha=0.8, width = 0.1)+
  stat_compare_means(paired = TRUE,tip.length = 0.02, size = 3, label.x.npc = "center")+
  scale_color_manual(values = c(blues9[4],blues9[8]))+
  theme_classic() + 
  stat_summary(fun.y= mean, geom="crossbar", width=0.2, color="black", alpha = 0.6)+
  stat_summary(geom = "linerange", fun.data = mean_sdl, fun.args = list(mult = 1),colour = "black", alpha = 0.6) +
  rremove("legend") + 
  xlab("")
ggsave(plot = plot, filename = "../Plots/SVM_dotplot_clonesvsrandom.svg", width = 3,height = 6)
write.csv(accuracy, "Accuracy_SVM_S1.csv", row.names = F)
#For S2
# ggsave(plot = plot, filename = "../Plots/SVM_dotplot_clonesvsrandom_S2.svg", width = 3,height = 6)
# write.csv(accuracy, "Accuracy_SVM_S2.csv", row.names = F)

#Plot combined dotplot
accuracy_S1 <- read.csv("Accuracy_SVM_S1.csv")
accuracy_S2 <- read.csv("Accuracy_SVM_S2.csv")
accuracy_comb <- rbind(data.frame(accuracy_S1,"SampleSet"= "S1" ), data.frame(accuracy_S2, "SampleSet"= "S2"))
plot <- ggplot(accuracy_comb, mapping = aes(x = Sample, y = Accuracy, color = SampleSet))+
  geom_jitter(alpha=0.8, width = 0.1)+
  stat_compare_means(aes(group = Sample), paired = TRUE,tip.length = 0.02, size = 3, label.x.npc = "center")+
  scale_color_manual(values = c(blues9[4],blues9[8]))+
  theme_classic() + 
  stat_summary(fun.y= mean, geom="crossbar", width=0.2, color="black", alpha = 0.6)+
  stat_summary(geom = "linerange", fun.data = mean_sdl, fun.args = list(mult = 1),colour = "black", alpha = 0.6) +
  xlab("")
ggsave(plot = plot, filename = "../Plots/SVM_dotplot_clonesvsrandom_combined.svg")
