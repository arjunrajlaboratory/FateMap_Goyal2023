########################################################################################
#Description: PCA for FateMap data
#Aim: PCA to visualize trajectory of tumor samples under treatment
#Author: Maalavika Pillai
#Version: 3
#Created on: 19/10/22
#Edit: Add labels, early to late line segment
########################################################################################
library(dplyr)
require(reshape2)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(DESeq2)
library(svglite)
library(ggrepel)

setwd("~/Library/CloudStorage/OneDrive-NorthwesternUniversity/Fatemap Revisions")
df <- read.delim("Data/Bulk Data/RC_RNASeq_JLPRRB_2_Normalized_meltedData.tsv", sep="\t") 

#Random matrix to have exactly 1 early sample
randMat <- function(df, n=2, n_early=1){
  rand_mat <- sample_n(df[df$group=="early",],n_early)
  rand_mat <- rbind(rand_mat, sample_n(df[df$group=="late",],n-n_early))
  return(rand_mat)
}

#Create dataframe with genes in rows and samples in columns
df_2 <- df %>% 
  reshape2::dcast(gene_name~sampleID, value.var = "counts", fun.aggregate = sum)
rownames(df_2) <- df_2$samppleID
df_2 <- df_2[,-1]
df_2 <- df_2[,paste0("RC-Sample",c(1:48))]

#Normalize data
samples <- read.delim("Data/Bulk Data/RC RNA seq - RNA conc..tsv", sep="\t")
samples$time <- "late"
samples$time[grepl("naive",samples$sample.name)] <- "Naive"
samples$time[grepl("early",samples$sample.name)] <- "early"
samples$Cell <- substr(samples$sample.name,1,5)
samples$Cell <- gsub("early.*","", samples$sample.name)
samples$Cell <- gsub("late.*","", samples$Cell)
samples$Cell <- gsub("naive.*","", samples$Cell)

#Set number of variable genes
var_genes = 200

#Remove naive samples
df_2 <- df_2[,!samples$time=="Naive"]
samples <- samples[!samples$time=="Naive",]

#Create DESeq object and normalize
DS_object <- DESeqDataSetFromMatrix(df_2,samples,~1)
DS_object <- DS_object[rowSums(counts(DS_object)) >= 10,]
vsd <- vst(DS_object,blind=FALSE)
vsd$time <- factor(vsd$time, levels = c("early","late"))

PCA <- plotPCA(vsd, intgroup=c("time"), returnData=FALSE, n = var_genes)+
  theme_classic()+
  geom_line(aes(group = vsd$Cell), color = "darkgrey")+
  #rremove("legend")+
  border()

PCA <- data.frame(PCA$data,'Cell'= samples$Cell)
PCA<- with(PCA, PCA[order(group),])

#set seed for reproducibility
set.seed(100)
seed <- sample(1:10000,length(unique(PCA$Cell)))

dist_mat <- data.frame()
for (i in unique(PCA$Cell)) {
  pair_mat <- PCA[PCA$Cell==i,1:3]
  n <- nrow(pair_mat)
  if(n>1&all(c("early","late") %in% pair_mat$group)){
   #set.seed(seed[which(unique(PCA$Cell)==i)])
    rand_mat <- randMat(PCA[1:3],n)
    dist_mat <- rbind(dist_mat,data.frame('Distance' = dist(pair_mat[,-3])[1:n-1],'Type'= rep("Paired",n-1), 'Index'= i))
    dist_mat <- rbind(dist_mat,data.frame('Distance' = dist(rand_mat[,-3])[1:n-1],'Type'=  rep("Random",n-1), 'Index'= i))
  }
}

plot <- ggplot(dist_mat, aes(x=Type, y=Distance))+ 
  geom_jitter(alpha=0.3)+
  stat_compare_means(paired = TRUE, label.x.npc = "center",tip.length = 0.02, size = 5)+
  scale_color_manual(values = c("#A00000","darkgrey" ))+
  theme_classic((base_size = 22)) + 
  geom_boxplot(alpha=0, outlier.shape = NA, position = position_dodge(1), width = 0.8)+
 # geom_line(aes(group=Index), alpha =0.1)+
  rremove("legend") + 
  ylim(0,90)+
  xlab("")
ggsave(file=paste0("Plots/Bulk/revisionRCDrift_dist2PCs_", var_genes, "var.svg"), plot=plot, width=6, height=8)

#Manually extract PCA
df <- assay(vsd)

#For top 500 variable genes
sel = order(apply(df, 1, var), decreasing=TRUE)[1:var_genes]
df <- df[sel,]

#Calculate PCs
PCA <- prcomp(t(df))
axes=1:2
eig <- get_eigenvalue(PCA)[axes,2]
eig <- round(eig,2)
variance <- data.frame(PC= paste0(1:length(PCA$sdev)),
                       var_explained=get_eigenvalue(PCA)[,2])
variance$PC <- factor(variance$PC, levels = paste0(1:length(PCA$sdev)))

#Plotting scree plot
plot <- variance %>%
  ggplot(aes(x=PC,y=cumsum(var_explained),group=1))+
  geom_point(size=2)+
  geom_line()+
  theme_classic()+
  labs(title="Scree plot: PCA on scaled data")+
  xlab("Principle component")+
  ylab("Cumulative Percentage variance explained")+
  border()
ggsave(file=paste0("Plots/Bulk/revisionRCDrift_screePlot_", var_genes, "var.svg"), plot=plot, width=10, height=8)
#Extract values for 25PCs
ind <- facto_summarize(PCA, element = "ind", result = "coord", axes = 1:25)
PCA <- data.frame(ind[,c(paste0("Dim.",1:25))], 'Cell'= samples$Cell, group = samples$time)
PCA<- with(PCA, PCA[order(group),])
segment <- data.frame()
for (i in unique(PCA$Cell)) {
  PCA_sub <- PCA[PCA$Cell==i,]
  if(all(c("early","late") %in% PCA_sub$group)){
    segment_sub <-PCA_sub[PCA_sub$group=="late",c("Dim.1","Dim.2")]
    segment_sub$Dim.1.early <-  PCA_sub[PCA_sub$group=="early",c("Dim.1")]
    segment_sub$Dim.2.early <-  PCA_sub[PCA_sub$group=="early",c("Dim.2")]
    segment <- rbind(segment, segment_sub)
  }
} 

#Plotting PCA - sanity check to make sure it looks same as one from PlotPCA
PCA_main <- PCA 
PCA$Cell[PCA$group=="late"] <- ""
plot <- ggplot(PCA,aes(x=Dim.1, y=Dim.2))+
  geom_text_repel(data= PCA, aes(x=Dim.1, y=Dim.2, label = Cell))+
  geom_point( data= PCA,color=factor(PCA$group, levels = c("early","late"), labels = c(blues9[4],blues9[8])), size=3)+
  geom_segment(data = segment,aes(x = segment$Dim.1.early, xend = segment$Dim.1 , y = segment$Dim.2.early, yend = segment$Dim.2 ),color = "darkgrey")+
  theme_classic()+
  xlab(paste0("PC1 (",eig[1],"%)"))+
  ylab(paste0("PC2 (",eig[2],"%)"))+
  border()
ggsave(file=paste0("Plots/Bulk/revisionRCDrift_PCA_", var_genes, "var.svg"), plot=plot, width=10, height=8)

plot <- ggplot(PCA,aes(x=Dim.1, y=Dim.2))+
  geom_point( data= PCA,color=factor(PCA$group, levels = c("early","late"), labels = c(blues9[4],blues9[8])), size=3)+
  geom_segment(data = segment,aes(x = segment$Dim.1.early, xend = segment$Dim.1 , y = segment$Dim.2.early, yend = segment$Dim.2 ),color = "darkgrey")+
  theme_classic()+
  xlab(paste0("PC1 (",eig[1],"%)"))+
  ylab(paste0("PC2 (",eig[2],"%)"))+
  border()
ggsave(file=paste0("Plots/Bulk/revisionRCDrift_PCA_noLabels_", var_genes, "var.svg"), plot=plot, width=10, height=8)
PCA <- PCA_main

#set seed for reproducibility
set.seed(12)
seed <- sample(1:10000,length(unique(PCA$Cell)))

#Calculate distance based on 25 PCs
dist_mat <- data_frame()
for (i in unique(PCA$Cell)) {
  pair_mat <- as.data.frame(PCA[PCA$Cell==i,])
  n <- nrow(pair_mat)
  if(n>1&all(c("early","late") %in% pair_mat$group)){
    set.seed(seed[which(unique(PCA$Cell)==i)])
    rand_mat <- as.data.frame(randMat(PCA,n))
    dist_mat <- rbind(dist_mat,data.frame('Distance' = dist(pair_mat[,paste0("Dim.",1:25)])[1:n-1],'Type'= rep("Paired",n-1), 'Index'= i))
    dist_mat <- rbind(dist_mat,data.frame('Distance' = dist(rand_mat[,paste0("Dim.",1:25)])[1:n-1],'Type'=  rep("Random",n-1), 'Index'= i))
  }
}

#Plotting barplot
plot <- ggplot(dist_mat, aes(x=Type, y=Distance))+ 
  geom_jitter(alpha=0.3)+
  stat_compare_means(paired = TRUE, label.x.npc = "center",tip.length = 0.02, size = 5)+
  scale_color_manual(values = c("#A00000","darkgrey" ))+
  theme_classic((base_size = 22)) + 
  geom_boxplot(alpha=0, outlier.shape = NA, position = position_dodge(1), width = 0.8)+
  rremove("legend") + 
  ylim(0,90)
  xlab("")
ggsave(file=paste0("Plots/Bulk/revisionRCDrift_dist25PCs_", var_genes, "var.svg"), plot=plot, width=6, height=8)


