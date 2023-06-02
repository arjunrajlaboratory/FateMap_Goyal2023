#####################NOTE#####################
#Cleaned up version of V2: /Users/yogesh/Dropbox (RajLab)/FateMap/paper/scripts_sept2021_onward/PCA_comparison.R
##############################################
########################################################################################
#Description: PCA variance analysis FM06 
#Author: Yogesh Goyal
#Edited by: Maalavika Pillai
#Version: 2.2
#Edited on: 10/24/22
#Funtion created to generate plot, includees analysis for pre-resistant samples
########################################################################################

library(plyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(svglite)
library(Seurat)
library(dplyr)

# Load the relevant datasets (filteredbarcode matrix files)

home1Directory <-"~/Library/CloudStorage/OneDrive-NorthwesternUniversity/Fatemap Revisions/Data/Drug treated and naive/FM01/"
home2Directory <-"~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2021_FM06/extracted/FM06/"
home3Directory <-"~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2021_FM06/extracted/PrimedCells/"
plotDirectory <-"~/Dropbox (RajLab)/FateMap/REVISIONS/FateMapAnalysis_Jonas/2021_FM06/plot/"

PCAvarPlot <- function(FM, shuffles = 100, label = "FM"){
  require(plyr)
  require(ggplot2)
  require(dplyr)
  #follow zeehio's response on https://github.com/satijalab/seurat/issues/982/
  FM <- RunPCA(FM,assay='integrated',slot="scale.data")
  mat_fm <- GetAssayData(FM,assay='integrated',slot="scale.data") #rows are genes, columns are cells
  pca_fm<- FM[['pca']]
  total_variance_fm <- sum(matrixStats::rowVars(mat_fm)) #find variance estimate per row, i.e. per gene
  
  eigValues_fm=(pca_fm@stdev)^2
  varExplained_fm = as_tibble(eigValues_fm/total_variance_fm)%>% #sum is 0.28 (<<1 bc this is only first 50 PCs maybe??)
    transmute(varExplained = value) %>%
    mutate(PC = row_number(),cumulative = cumsum(varExplained))
  
  fmvar = ggplot(varExplained_fm,aes(x=PC))+
    geom_bar(aes(y=varExplained),stat='identity')+
    xlab("PC")+
    ylab("fraction of variance explained")+
    ggtitle("Resistant sample variance explained")+
    theme_bw()
  print(fmvar)
  
  randomVarBootStrap <- data.frame()
  set.seed(100)
  seed <- sample(1:(shuffles*10), shuffles)
  for (bootstrap in 1:shuffles) {
    mat_fm_random <- mat_fm 
    set.seed(seed[bootstrap])
    for (index_row in 1:nrow(mat_fm_random)){
      mat_fm_random[index_row,] <- sample(mat_fm_random[index_row,], ncol(mat_fm))
    }
    fm_random <- SetAssayData(
      object = FM,
      slot = "scale.data",
      new.data = mat_fm_random,
      assay = "integrated"
    )
    
    fm_random <- RunPCA(fm_random, assay='integrated',features = VariableFeatures(object = fm_random), verbose = FALSE)
    
    pca_fm_random <- fm_random[['pca']]
    total_variance_fm_random <- sum(matrixStats::rowVars(mat_fm_random)) #find variance estimate per row, i.e. per gene
    
    eigValues_fm_random=(pca_fm_random@stdev)^2
    varExplained_fm_random = as_tibble(eigValues_fm_random/total_variance_fm_random)%>% #sum is 0.729 (<1 bc this is only first 50 PCs maybe)
      transmute(varExplained = value) %>%
      mutate(PC = row_number(),cumulative = cumsum(varExplained), bootStrap=bootstrap)
    randomVarBootStrap <- rbind(randomVarBootStrap , varExplained_fm_random)
  }
  
  data_summary <- function(data, varname, groupnames){
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- plyr::rename(data_sum, c("mean" = varname))
    return(data_sum)
  }
  names(randomVarBootStrap)
  cutPoint <- vector()
   for (i in 1:shuffles ) {
    randomVar_sub <- randomVarBootStrap[randomVarBootStrap$bootStrap == i,]
    PC <- randomVar_sub$PC[which(randomVar_sub$varExplained > varExplained_fm$varExplained )]
    print(PC)
    cutPoint[i] <- min(PC)
  }
  randomVarBootStrap <- data_summary(data = randomVarBootStrap, groupnames = c("PC"), varname = "varExplained")
  
  
  fm_varrandom = ggplot(randomVarBootStrap, aes(x=PC, y=varExplained)) + 
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=varExplained-sd, ymax=varExplained+sd), width=.2,
                  position=position_dodge(.9)) +
    xlab("PC")+
    ylab("fraction of variance explained")+
    ggtitle(paste0(label," randomized variance explained"))+
    theme_bw()
  print(fm_varrandom)
  
  #add to original iPS plot
  fmvar_line = fmvar +
    geom_line(data=randomVarBootStrap, aes(x=PC, y=varExplained),stat='identity',color='#F8766D')+
    geom_errorbar(data= randomVarBootStrap, aes(ymin=varExplained-sd, ymax=varExplained+sd), width=.2,
                  position=position_dodge(.9),color='#F8766D') +
    xlim(0,50) +
    theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1))
  
  return(list(fmvar_line, fmvar, randomVarBootStrap ,varExplained_fm , cutPoint ))
}


fm06_sample1_2 = readRDS(file = paste0(home2Directory,'s1s2Naive_ScTransform_50pcs_filter.rds'))  ##to read if needed.
f1 <- PCAvarPlot(fm06_sample1_2 , shuffles = 100)  
write.table(f1[[3]], paste0(home2Directory,"RandomizedPCAVariance.tsv"), sep = "\t", quote = F, row.names = F)
write.table(f1[[2]]$data, paste0(home2Directory,"OriginalPCAVariance.tsv"), sep = "\t", quote = F, row.names = F)
ggsave(plot = f1[[1]],  filename = paste0(plotDirectory, 'fm06_resistant_PCA_100bootstraps.svg'), width = 8, height = 7)

fm01_sample1_2 = readRDS(file = paste0(home1Directory,'s1s2ScTransform_50pcs_filter.rds')) 
f2 <- PCAvarPlot(fm01_sample1_2, shuffles = 100)  
write.table(f2[[3]], paste0(home1Directory,"RandomizedPCAVariance.tsv"), sep = "\t", quote = F, row.names = F)
write.table(f2[[2]]$data, paste0(home1Directory,"OriginalPCAVariance.tsv"), sep = "\t", quote = F, row.names = F)
ggsave(plot = f2[[1]],  filename = paste0(plotDirectory, 'fm01_resistant_PCA_100bootstraps.svg'), width = 8, height = 7)

fm06_subsamp_0.4 = readRDS(file = paste0(home3Directory,'FM06_preresistant_0.4.rds'))  ##to read if needed.
f3 <- PCAvarPlot(fm06_subsamp_0.4 , shuffles = 100)  
write.table(f3[[3]], paste0(home3Directory,"RandomizedPCAVariance.tsv"), sep = "\t", quote = F, row.names = F)
write.table(f3[[2]]$data, paste0(home3Directory,"OriginalPCAVariance.tsv"), sep = "\t", quote = F, row.names = F)
ggsave(plot = f3[[1]],  filename = paste0(plotDirectory, 'fm06_pre-resistant_0.04_PCA_100bootstraps.svg'), width = 8, height = 7)

# fm06_subsamp_0.5 = readRDS(file = paste0(home3Directory,'FM06_preresistant_0.5.rds'))  ##to read if needed.
# f4 <- PCAvarPlot(fm06_subsamp_0.5 , shuffles = 100)  
remove(fm01_sample1_2, fm06_subsamp_0.4, fm06_sample1_2)
save.image( paste0(home3Directory,"All3_images.RData"))


#Combine data from f1,f2,f3
fm06_var <- data.frame(rbind(data.frame(f1[[4]][1:2],'sd' = NA ,'Type' = "original"),data.frame(f1[[3]][1:3], 'Type' = "random")), "Sample" = "Naive")
fm01_var <- data.frame(rbind(data.frame(f2[[4]][1:2], 'sd' = NA ,'Type' = "original"),data.frame(f2[[3]][1:3], 'Type' = "random")), "Sample" = "Resistant")
primedCell_var <- data.frame(rbind(data.frame(f3[[4]][1:2], 'sd' =NA ,'Type' = "original"),data.frame(f3[[3]][1:3], 'Type' = "random")),"Sample" = "PrimedCell")
var_all <- rbind(rbind(fm06_var, fm01_var), primedCell_var)

#Plot all together as intersecting line plots
plot <- ggplot(data = var_all) +
  geom_point(aes(x = PC, y = varExplained, color = Type))+
  geom_errorbar(data = var_all,aes(x = PC, ymin = varExplained - sd, ymax = varExplained +sd))+
  facet_wrap(.~Sample, ncol =1, scales = "free")+
  theme_classic() 
ggsave(filename = paste0(plotDirectory, "PCAvarianceRandom.svg"), plot= plot, height = 10, width = 8)


#Plot cut-off ofr each population
cutOff_stat <- rbind(rbind(data.frame("Cut" = f1[[5]], "Sample" = "Naive"), data.frame("Cut" = f2[[5]], "Sample" = "Resistant")),data.frame("Cut" = f3[[5]], "Sample" = "Primed Cell"))
stat <- compare_means(Cut ~ Sample, data = cutOff_stat )
stat <- stat%>%
  mutate(y.position = c(45, 35, 48))
cutOff <- data.frame("PC" = c(mean(f1[[5]]), mean(f2[[5]]), mean(f3[[5]])), "sd" = c(sd(f1[[5]]), sd(f2[[5]]), sd(f3[[5]])), "Sample" = c("Naive", "Resistant", "Primed Cell"))
cutOff$Sample <- factor(cutOff$Sample, levels = c("Resistant", "Naive", "Primed Cell"))
plot <- ggplot(cutOff, aes(x=factor(Sample), y=PC, group = 1)) + 
  geom_point()+
  geom_line(color = "black")+
  geom_errorbar(aes(ymin=PC-sd, ymax=PC+sd), width=.2,
                position=position_dodge(0.05))+
  stat_pvalue_manual(stat)+
  NoLegend()+
  theme_classic()+
  xlab("")+
  ylab("Principal Components needed to explain variance")
ggsave(filename = paste0(plotDirectory, "PCAvarianceRandom_points.svg"), plot= plot, height = 6, width = 6) 
