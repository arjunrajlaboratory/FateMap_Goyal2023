#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/gDNA_BCSeq/20210901_Mouse_Run2/analyzed/Expt1/scripts/2021.09.16_Mouse01_barcodeAnalysis.R
##############################################

library(tidyverse)
library(reshape2)
library(ggridges)

homeDirectory = '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/gDNA_BarcodeData/mouseExpts/'
plotDirectory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/plotData/analysisPlots/'
makeOverlapTable <- function(a,b){
  rankList <- c(50, 100, 250, 500, 1000)
  
  overlapTable <- data.frame(matrix(NA, nrow = 5, ncol = 5))
  colnames(overlapTable) <-
    c("A50", "A100", "A250", "A500", "A1000")
  rownames(overlapTable) <-
    c("B50", "B100", "B250", "B500", "B1000")
  
  for (i in rankList) {
    for (j in rankList) {
      aTopTable <- top_n(a, i, V2)
      bTopTable <- top_n(b, j, V2)
      
      overlapTemp <- inner_join(aTopTable, bTopTable, by = "V1")
      overlapVal <- nrow(overlapTemp) / min(nrow(aTopTable), nrow(bTopTable))
      overlapTable[paste("B", i, sep = ""), paste("A", j, sep = "")] = round(overlapVal,3)
    }
  }
  
  dfm <- melt(as.matrix(overlapTable))
  
  table <- ggplot(dfm, aes(
    x = Var2,
    y = Var1,
    label = value,
    fill = value
  )) +
    geom_tile(show.legend = FALSE) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red4", limits = c(0, 1), midpoint = 0.5) +
    geom_text(color = "black") +
    #xlab("Top n Barcodes in A") +
    #ylab("Top n Barcodes in B") + 
    theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))
  
  return(table)
}
cutoff = c(600, 800)

summary = matrix(, nrow = 3*length(cutoff), ncol = 7)

sampleDir <- paste0(homeDirectory, 'Expt1/')
sampleFolders <- list.files(path=sampleDir,pattern='*d8.txt',recursive=TRUE)

sampleTables <- list(rep("NA",length(sampleFolders)))
for (i in 1:length(sampleFolders)) {
  sampleTable = read.table(paste0(sampleDir, sampleFolders[i]), header = FALSE, sep = "\t")
  sampleTables[[i]] <- sampleTable
}


###############################################################################################################################
###########################################DRUG DEPENDENT vs INDEPENDENT LINEAGES##############################################
cond2 = "305 (2)"
cond1 = "301 (1)"
filterThreshold <- 0
x = 1; 
y = 2;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*(10^6))
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*(10^6))

overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x, V2.y=V2.y) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
foldChangecutoff = 2.5

for (i in 1:length(cutoff)) {
drugDepX <- overlapPLX %>% filter(foldchange < -foldChangecutoff | foldchange > foldChangecutoff) %>% filter(V2.x>cutoff[i] & V2.y<cutoff[i])
drugDepY <- overlapPLX %>% filter(foldchange < -foldChangecutoff | foldchange > foldChangecutoff) %>% filter(V2.y>cutoff[i] & V2.x<cutoff[i])
drugInd <- overlapPLX %>% filter(foldchange > -foldChangecutoff & foldchange < foldChangecutoff) %>% filter(V2.x>cutoff[i] & V2.y>cutoff[i])

plot <- ggplot() +
  geom_point(aes(x=overlapPLX$V2.x,y=overlapPLX$V2.y),alpha=0.4, size = 4, shape = 16) +
  geom_point(aes(x=drugInd$V2.x,y=drugInd$V2.y),col="orange", size = 4, shape = 16) +
  geom_point(aes(x=drugDepY$V2.x,y=drugDepY$V2.y),col="springgreen3", size = 4, shape = 16) +
  geom_point(aes(x=drugDepX$V2.x,y=drugDepX$V2.y),col="turquoise3", size = 4, shape = 16) +
  geom_abline(slope=1,intercept=0) +
  xlab(paste0(cond1)) + ylab(paste0(cond2)) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 2), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot = plot, file = paste0(plotDirectory, 'BCMouse_independence',cond1,'Vs', cond2,'cutoff_',cutoff[i],'.svg'), width = 8, height = 7)

summary[i,] = c(nrow(drugDepX),nrow(drugDepY), nrow(drugInd), 100*(nrow(drugInd))/(nrow(drugInd) + nrow(drugDepY) + nrow(drugDepX)), (nrow(drugInd) + nrow(drugDepY) + nrow(drugDepX)), cutoff[i], "Expt1")
}

######
sampleDir <- paste0(homeDirectory, 'Expt2/')
sampleFolders <- list.files(path=sampleDir,pattern='*d8.txt',recursive=TRUE)

sampleTables <- list(rep("NA",length(sampleFolders)))
for (i in 1:length(sampleFolders)) {
  sampleTable = read.table(paste0(sampleDir, sampleFolders[i]), header = FALSE, sep = "\t")
  sampleTables[[i]] <- sampleTable
}

###############################################################################################################################
###########################################DRUG DEPENDENT vs INDEPENDENT LINEAGES##############################################
cond2 = "306 (2)"
cond1 = "310 (1)"
filterThreshold <- 0
x = 1; 
y = 2;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*(10^6))
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*(10^6))

overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x, V2.y=V2.y) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
foldChangecutoff = 2.5

for (i in 1:length(cutoff)) {
  drugDepX <- overlapPLX %>% filter(foldchange < -foldChangecutoff | foldchange > foldChangecutoff) %>% filter(V2.x>cutoff[i] & V2.y<cutoff[i])
  drugDepY <- overlapPLX %>% filter(foldchange < -foldChangecutoff | foldchange > foldChangecutoff) %>% filter(V2.y>cutoff[i] & V2.x<cutoff[i])
  drugInd <- overlapPLX %>% filter(foldchange > -foldChangecutoff & foldchange < foldChangecutoff) %>% filter(V2.x>cutoff[i] & V2.y>cutoff[i])

  summary[i+1*length(cutoff),] = c(nrow(drugDepX),nrow(drugDepY), nrow(drugInd), 100*(nrow(drugInd))/(nrow(drugInd) + nrow(drugDepY) + nrow(drugDepX)), (nrow(drugInd) + nrow(drugDepY) + nrow(drugDepX)), cutoff[i], "Expt2")
}

######
sampleDir <- paste0(homeDirectory, 'Expt5/')
sampleFolders <- list.files(path=sampleDir,pattern='*d8.txt',recursive=TRUE)

sampleTables <- list(rep("NA",length(sampleFolders)))
for (i in 1:length(sampleFolders)) {
  sampleTable = read.table(paste0(sampleDir, sampleFolders[i]), header = FALSE, sep = "\t")
  sampleTables[[i]] <- sampleTable
}


###############################################################################################################################
###########################################DRUG DEPENDENT vs INDEPENDENT LINEAGES##############################################
cond2 = "324 (2)"
cond1 = "325 (1)"
filterThreshold <- 0
x = 1; 
y = 2;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*(10^6))
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*(10^6))

overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x, V2.y=V2.y) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
foldChangecutoff = 2.5

for (i in 1:length(cutoff)) {
  drugDepX <- overlapPLX %>% filter(foldchange < -foldChangecutoff | foldchange > foldChangecutoff) %>% filter(V2.x>cutoff[i] & V2.y<cutoff[i])
  drugDepY <- overlapPLX %>% filter(foldchange < -foldChangecutoff | foldchange > foldChangecutoff) %>% filter(V2.y>cutoff[i] & V2.x<cutoff[i])
  drugInd <- overlapPLX %>% filter(foldchange > -foldChangecutoff & foldchange < foldChangecutoff) %>% filter(V2.x>cutoff[i] & V2.y>cutoff[i])
  
  summary[i+2*length(cutoff),] = c(nrow(drugDepX),nrow(drugDepY), nrow(drugInd), 100*(nrow(drugInd))/(nrow(drugInd) + nrow(drugDepY) + nrow(drugDepX)), (nrow(drugInd) + nrow(drugDepY) + nrow(drugDepX)), cutoff[i], "Expt5")
}

write.table(summary, file=paste0(plotDirectory,'BCMouse_summaryTable.tsv'), col.names = TRUE, sep='\t')


