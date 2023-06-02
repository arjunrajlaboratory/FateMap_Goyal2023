#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/gDNA_BCSeq/20201112_gDNA_BCSeq/20201112_BC03/scripts/2020.11.12_BC03_barcodeAnalysis.R
##############################################

library(tidyverse)
library(reshape2)
library(ggridges)

sampleDir <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/gDNA_BarcodeData/BC03/analyzed/'
sampleFolders <- list.files(path=sampleDir,pattern='*d8.txt',recursive=TRUE)
plot3Directory <- "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/plotData/analysisPlots/"

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

#[1] "BC03_1_1uM_PLX_1/1_1uM_PLX_1_clusteredBarcodeUMIs_d8.txt" 
#[2] "BC03_1_1uM_PLX_2/1_1uM_PLX_2_clusteredBarcodeUMIs_d8.txt"            
#[3] "BC03_2_100nM_PLX/2_100nM_PLX_clusteredBarcodeUMIs_d8.txt"                 
#[4] "BC03_3_100nM_PLX/3_100nM_PLX_clusteredBarcodeUMIs_d8.txt"
#[5] "BC03_4_1uM_PLX/4_1uM_PLX_clusteredBarcodeUMIs_d8.txt"
#[6] "BC03_Naive/Naive_clusteredBarcodeUMIs_d8.txt"

sampleTables <- list(rep("NA",length(sampleFolders)))
for (i in 1:length(sampleFolders)) {
  sampleTable = read.table(paste0(sampleDir, sampleFolders[i]), header = FALSE, sep = "\t")
  sampleTables[[i]] <- sampleTable
}

####Normalizing using spike-in barcodes
spikeBC1 = "TCCAGGTCCTCCTACTTGTACAACACCTTGTACAGCTGCTAGTGGTAGAAGAGGTACAACAACAACACGAGCATCATGAGGATCTACAGCATCAAGAACA"
spikeBC2 = "ACGTTGTGCATGACCTTGATCACCAGCTCGATGTCGAACATCACGAGCTCGTTCTGCATCTGCAAGAACACCTCGTCCTTGAACTGCTCGACGTCCATGA"  
nBC1tot = c(25000,25000,50000, 50000,25000,10000) #1-4
nBC2tot = c(5000,5000,10000,10000,5000,2000)    #1-4
nFraction = c(0.41, 0.41, 0.33, 0.31, 0.52, 0.51)    #1-4
nBC1 = nBC1tot*nFraction
nBC2 = nBC2tot*nFraction

standardTableAll = list()
lmr = list()

for (i in 1:length(sampleFolders)) {
  standardTable = sampleTables[[i]] %>% filter(sampleTables[[i]]$V1 == spikeBC1 | sampleTables[[i]]$V1 == spikeBC2) %>% mutate(sampleNum = i, nBC = c(nBC1[i],nBC2[i]))
  standardTable = standardTable %>% mutate(ratio = standardTable$V2[1]/standardTable$V2[2])
  lmr[i] = coef(lm(standardTable$nBC ~ standardTable$V2 - 1))
  if(is.null(dim(standardTableAll))){
    standardTableAll = standardTable
  } else {
    standardTableAll = bind_rows(standardTableAll, standardTable)
  }
}
#standardTableAll =  as_tibble(standardTableAll %>% add_row(V1 = 'noBC', V2 = 0, sampleNum = 0, nBC = 0, ratio = NaN))

for (i in 1:length(sampleFolders)) {
  sampleTables[[i]] = sampleTables[[i]] %>% dplyr::filter(sampleTables[[i]]$V1 != spikeBC1)
  sampleTables[[i]] = sampleTables[[i]] %>% dplyr::filter(sampleTables[[i]]$V1 != spikeBC2)
}


###############################################################################################################################
###########################################DRUG DEPENDENT vs INDEPENDENT LINEAGES##############################################
###########################################OverCutoff##############################################

# (2,3) 100nM PLX vs 1uM PLX
cond2 = "100nM PLX (3)"
cond1 = "1uM PLX (2)"
filterThreshold <- 0
x = 2; 
y = 3;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x, V2.y=V2.y) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)

cutoff = c(100,200,300,400,500,600,700,800,900,1000)
foldChangecutoffDep = 2.5
foldChangecutoffInd = 1.5

ratioDepXDepY = c()
ratioDepXInd = c()
ratioDepYInd = c()
ratioIndDepYIndDepX = c()

for (i in 1:length(cutoff)) {
  drugDepX <- overlapPLX %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(V2.x>cutoff[i] & V2.y<cutoff[i])
  drugDepY <- overlapPLX %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(V2.y>cutoff[i] & V2.x<cutoff[i])
  drugInd <- overlapPLX %>% filter(foldchange > -foldChangecutoffInd & foldchange < foldChangecutoffInd) %>% filter(V2.x>cutoff[i] & V2.y>cutoff[i])
  
  ratioDepXDepY[i] = nrow(drugDepY)/nrow(drugDepX)
  ratioDepXInd[i] = nrow(drugInd)/nrow(drugDepX)
  ratioDepYInd[i] = nrow(drugInd)/nrow(drugDepY)
  ratioIndDepYIndDepX[i] = (nrow(drugInd) + nrow(drugDepY))/(nrow(drugInd) + nrow(drugDepX))
  
  plot <- ggplot() +
    geom_point(aes(x=overlapPLX$V2.x,y=overlapPLX$V2.y),alpha=0.4, size = 4, shape = 16) +
    geom_point(aes(x=drugInd$V2.x,y=drugInd$V2.y),col="orange", size = 4, shape = 16) +
    geom_point(aes(x=drugDepY$V2.x,y=drugDepY$V2.y),col="springgreen3", size = 4, shape = 16) +
    geom_point(aes(x=drugDepX$V2.x,y=drugDepX$V2.y),col="turquoise3", size = 4, shape = 16) +
    geom_abline(slope=1,intercept=0) +
    xlab(paste0(cond1)) + ylab(paste0(cond2)) +
    theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 2), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  ggsave(plot = plot, file = paste0(plot3Directory, 'BC03_independence',cond1,'Vs', cond2,'cutoff_',cutoff[i],'.svg'), width = 8, height = 7)
}
#########################################################################################

#########################################################################################
##############Analysis to see where the green cells go; for Supplementary ###################################
#########################################################################################
cond2 = "100nM PLX (3)"
cond1 = "1uM PLX (2)"
filterThreshold <- 0
x = 2; 
y = 3;
cutoff = 400;

a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x+1, V2.y=V2.y+1) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
drugDepY_1 <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>cutoff & V2.x<cutoff) %>% select(V1) %>% unique()

cond2 = "100nM PLX (4)"
cond1 = "1uM PLX (5)"
filterThreshold <- 0
x = 5; 
y = 4;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x+1, V2.y=V2.y+1) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
drugDepY_2 <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>cutoff & V2.x<cutoff) %>% select(V1) %>% unique()

cond2 = "100nM PLX (3)"
cond1 = "1uM PLX (5)"
filterThreshold <- 0
x = 5; 
y = 3;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x+1, V2.y=V2.y+1) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
drugDepY_3 <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>cutoff & V2.x<cutoff) %>% select(V1) %>% unique()

cond2 = "100nM PLX (4)"
cond1 = "1uM PLX (2)"
filterThreshold <- 0
x = 2; 
y = 4;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x+1, V2.y=V2.y+1) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
drugDepY_4 <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>cutoff & V2.x<cutoff) %>% select(V1) %>% unique()

# drugDepY_12 = inner_join(drugDepY_1,drugDepY_2)
# drugDepY_34 = inner_join(drugDepY_3,drugDepY_4)

drugDepY_All <- bind_rows(drugDepY_1,drugDepY_2,drugDepY_3,drugDepY_4) %>% select(V1) %>% unique()
# drugDepY_AllB <- bind_rows(drugDepY_12,drugDepY_34) %>% unique()

cond2 = "100nM PLX (1) only-Dependent"
cond1 = "100nM PLX (4) only-Dependent"
filterThreshold <- 0
x = 3; 
y = 4;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x+1, V2.y=V2.y+1) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)

drugDepX <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.x>cutoff & V2.y<cutoff)
drugDepY <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>cutoff & V2.x<cutoff)
drugInd <- overlapPLX %>% filter(foldchange > -1.5 & foldchange < 1.5) %>% filter(V2.x>cutoff & V2.y>cutoff)
drugDepY_All_joined = inner_join(overlapPLX,drugDepY_All, by = "V1")
# drugDepY_All_joinedB = inner_join(overlapPLX,drugDepY_AllB, by = "V1")

drugDepX1 <- drugDepY_All_joined %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.x>cutoff & V2.y<cutoff)
drugDepY1 <- drugDepY_All_joined %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>cutoff & V2.x<cutoff)
drugInd1 <- drugDepY_All_joined %>% filter(foldchange > -1.5 & foldchange < 1.5) %>% filter(V2.x>cutoff & V2.y>cutoff)

#####percentages for Venn plot:
percentDepX1 = nrow(drugDepX1)/(nrow(drugDepX1) + nrow(drugDepY1) + nrow(drugInd1))
percentDepY1 = nrow(drugDepY1)/(nrow(drugDepX1) + nrow(drugDepY1) + nrow(drugInd1))
percentdrugInd1 = nrow(drugInd1)/(nrow(drugDepX1) + nrow(drugDepY1) + nrow(drugInd1))

#####for Figure 3
plot <- ggplot() +
  geom_point(aes(x=overlapPLX$V2.x,y=overlapPLX$V2.y),alpha=0.4, size = 4, shape = 16) +
  geom_point(aes(x=drugDepY_All_joined$V2.x,y=drugDepY_All_joined$V2.y),col="orange", size = 4, shape = 16) +
  geom_point(aes(x=drugDepY1$V2.x,y=drugDepY1$V2.y),col="springgreen3", size = 4, shape = 16) +
  geom_point(aes(x=drugDepX1$V2.x,y=drugDepX1$V2.y),col="turquoise3", size = 4, shape = 16) +
  geom_abline(slope=1,intercept=0) +
  xlab(paste0(cond1)) + ylab(paste0(cond2))+
 # annotate(geom = "text", x = max(overlapPLX$V2.y)/2, y = max(overlapPLX$V2.y)/2, label = expression(italic("  y=x")), hjust = "left") + 
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 2), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot = plot, file = paste0(plot3Directory, 'BC03_independence',cond1,'Vs', cond2,'FORFIGUREOnlyGreen.svg'), width = 8, height = 7)

##########################################################################################
##############################Control with Naive Cells ###################################
##########################################################################################
cond1 = "Naive  (6)"
cond2 = "100nM PLX  (3)"
x = 6; 
y = 3;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x, V2.y=V2.y) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)

OnlyGreenSplit3 = inner_join(drugDepY_3,drugDepY_1) %>% unique()
dependentBarcodesSplit3 = as_tibble(OnlyGreenSplit3$V1) %>% dplyr::rename(V1 = value) 
dependentBarcodesSplit3Plot = inner_join(dependentBarcodesSplit3, overlapPLX, by ="V1")

scatterPlot<- ggplot() +
  geom_point(aes(x=overlapPLX$V2.x,y=overlapPLX$V2.y),alpha=0.4, size = 4, shape = 16) +
  geom_point(aes(x=dependentBarcodesSplit3Plot$V2.x,y=dependentBarcodesSplit3Plot$V2.y),col="springgreen3",size = 4, shape = 16) +
  xlab(paste0(cond1)) + ylab(paste0(cond2)) +
  geom_abline(slope=1,intercept=0) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 2), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot = scatterPlot, file = paste0(plot3Directory, 'BC03_scatter',cond1,'Vs', cond2,'.svg'), width = 8, height = 7)

####
cond1 = "Naive  (6)"
cond2 = "100nM PLX  (4)"
x = 6; 
y = 4;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x, V2.y=V2.y) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)

OnlyGreenSplit4 = inner_join(drugDepY_2,drugDepY_4) %>% unique()
dependentBarcodesSplit4 = as_tibble(OnlyGreenSplit4$V1) %>% dplyr::rename(V1 = value) 
dependentBarcodesSplit4Plot = inner_join(dependentBarcodesSplit4, overlapPLX, by ="V1")

scatterPlot<- ggplot() +
  geom_point(aes(x=overlapPLX$V2.x,y=overlapPLX$V2.y),alpha=0.4, size = 4, shape = 16) +
  geom_point(aes(x=dependentBarcodesSplit4Plot$V2.x,y=dependentBarcodesSplit4Plot$V2.y),col="springgreen3",size = 4, shape = 16) +
  xlab(paste0(cond1)) + ylab(paste0(cond2)) +
  geom_abline(slope=1,intercept=0) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 2), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot = scatterPlot, file = paste0(plot3Directory, 'BC03_scatter',cond1,'Vs', cond2,'.svg'), width = 8, height = 7)

#########################################################################################
##############Unique number of barcodes per conditions###################################
#########################################################################################

#[1] "BC03_1_1uM_PLX_1/1_1uM_PLX_1_clusteredBarcodeUMIs_d8.txt" 
#[2] "BC03_1_1uM_PLX_2/1_1uM_PLX_2_clusteredBarcodeUMIs_d8.txt"            
#[3] "BC03_2_100nM_PLX/2_100nM_PLX_clusteredBarcodeUMIs_d8.txt"                 
#[4] "BC03_3_100nM_PLX/3_100nM_PLX_clusteredBarcodeUMIs_d8.txt"
#[5] "BC03_4_1uM_PLX/4_1uM_PLX_clusteredBarcodeUMIs_d8.txt"
#[6] "BC03_Naive/Naive_clusteredBarcodeUMIs_d8.txt"

BC03_1_1uM_PLX_1 = length(unique(sampleTables[[1]]$V1))/nFraction[1]
BC03_1_1uM_PLX_2 = length(unique(sampleTables[[2]]$V1))/nFraction[2]
BC03_2_100nM_PLX = length(unique(sampleTables[[3]]$V1))/nFraction[3]
BC03_3_100nM_PLX = length(unique(sampleTables[[4]]$V1))/nFraction[4]
BC03_4_1uM_PLX = length(unique(sampleTables[[5]]$V1))/nFraction[5]
BC03_Naive = length(unique(sampleTables[[5]]$V1))


###########################################DRUG DEPENDENT vs INDEPENDENT LINEAGES##############################################
###########################################OverCutoff##############################################
###This constitutes BarcodeData on Figure 4 and supplementary for cutoff

# (2,3) 100nM PLX vs 1uM PLX
cond2 = "100nM PLX (3)"
cond1 = "1uM PLX (2)"
filterThreshold <- 0
x = 2; 
y = 3;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x, V2.y=V2.y) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)

cutoff = c(100,200,300,400,500,600,700,800,900,1000)
foldChangecutoffDep = 2.5
foldChangecutoffInd = 1.5

ratioDepXDepY = c()
ratioDepXInd = c()
ratioDepYInd = c()
ratioIndDepYIndDepX = c()

for (i in 1:length(cutoff)) {
  drugDepX <- overlapPLX %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(V2.x>cutoff[i] & V2.y<cutoff[i])
  drugDepY <- overlapPLX %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(V2.y>cutoff[i] & V2.x<cutoff[i])
  drugInd <- overlapPLX %>% filter(foldchange > -foldChangecutoffInd & foldchange < foldChangecutoffInd) %>% filter(V2.x>cutoff[i] & V2.y>cutoff[i])
  
  ratioDepXDepY[i] = nrow(drugDepY)/nrow(drugDepX)
  ratioDepXInd[i] = nrow(drugInd)/nrow(drugDepX)
  ratioDepYInd[i] = nrow(drugInd)/nrow(drugDepY)
  ratioIndDepYIndDepX[i] = (nrow(drugInd) + nrow(drugDepY))/(nrow(drugInd) + nrow(drugDepX))
  
  plot <- ggplot() +
      geom_point(aes(x=overlapPLX$V2.x,y=overlapPLX$V2.y),alpha=0.4, size = 4, shape = 16) +
      geom_point(aes(x=drugInd$V2.x,y=drugInd$V2.y),col="orange", size = 4, shape = 16) +
      geom_point(aes(x=drugDepY$V2.x,y=drugDepY$V2.y),col="springgreen3", size = 4, shape = 16) +
      geom_point(aes(x=drugDepX$V2.x,y=drugDepX$V2.y),col="turquoise3", size = 4, shape = 16) +
      geom_abline(slope=1,intercept=0) +
      xlab(paste0(cond1)) + ylab(paste0(cond2)) +
      theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 2), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  ggsave(plot = plot, file = paste0(plot3Directory, 'BC03_independence',cond1,'Vs', cond2,'cutoff_',cutoff[i],'.svg'), width = 8, height = 7)
}
tableReplicate1 = tibble(cutoff = cutoff, ratioDepXDepY = ratioDepXDepY, ratioDepXInd = ratioDepXInd, ratioDepYInd = ratioDepYInd,ratioIndDepYIndDepX = ratioIndDepYIndDepX, rep = "rep1")

# (2,3) 100nM PLX vs 1uM PLX
cond2 = "100nM PLX (4)"
cond1 = "1uM PLX (5)"
filterThreshold <- 0
x = 5; 
y = 4;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x, V2.y=V2.y) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)

cutoff = c(100,200,300,400,500,600,700,800,900,1000)
foldChangecutoffDep = 2.5
foldChangecutoffInd = 1.5

ratioDepXDepY = c()
ratioDepXInd = c()
ratioDepYInd = c()
ratioIndDepYIndDepX = c()

for (i in 1:length(cutoff)) {
  drugDepX <- overlapPLX %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(V2.x>cutoff[i] & V2.y<cutoff[i])
  drugDepY <- overlapPLX %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(V2.y>cutoff[i] & V2.x<cutoff[i])
  drugInd <- overlapPLX %>% filter(foldchange > -foldChangecutoffInd & foldchange < foldChangecutoffInd) %>% filter(V2.x>cutoff[i] & V2.y>cutoff[i])
  
  ratioDepXDepY[i] = nrow(drugDepY)/nrow(drugDepX)
  ratioDepXInd[i] = nrow(drugInd)/nrow(drugDepX)
  ratioDepYInd[i] = nrow(drugInd)/nrow(drugDepY)
  ratioIndDepYIndDepX[i] = (nrow(drugInd) + nrow(drugDepY))/(nrow(drugInd) + nrow(drugDepX))
  
}
tableReplicate2 = tibble(cutoff = cutoff, ratioDepXDepY = ratioDepXDepY, ratioDepXInd = ratioDepXInd, ratioDepYInd = ratioDepYInd, ratioIndDepYIndDepX = ratioIndDepYIndDepX, rep = "rep2")

#########
cond2 = "100nM PLX (4)"
cond1 = "1uM PLX (2)"
filterThreshold <- 0
x = 2; 
y = 4;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x, V2.y=V2.y) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)

cutoff = c(100,200,300,400,500,600,700,800,900,1000)
foldChangecutoffDep = 2.5
foldChangecutoffInd = 1.5

ratioDepXDepY = c()
ratioDepXInd = c()
ratioDepYInd = c()
ratioIndDepYIndDepX = c()

for (i in 1:length(cutoff)) {
  drugDepX <- overlapPLX %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(V2.x>cutoff[i] & V2.y<cutoff[i])
  drugDepY <- overlapPLX %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(V2.y>cutoff[i] & V2.x<cutoff[i])
  drugInd <- overlapPLX %>% filter(foldchange > -foldChangecutoffInd & foldchange < foldChangecutoffInd) %>% filter(V2.x>cutoff[i] & V2.y>cutoff[i])
  
  ratioDepXDepY[i] = nrow(drugDepY)/nrow(drugDepX)
  ratioDepXInd[i] = nrow(drugInd)/nrow(drugDepX)
  ratioDepYInd[i] = nrow(drugInd)/nrow(drugDepY)
  ratioIndDepYIndDepX[i] = (nrow(drugInd) + nrow(drugDepY))/(nrow(drugInd) + nrow(drugDepX))
  
}
tableReplicate3 = tibble(cutoff = cutoff, ratioDepXDepY = ratioDepXDepY, ratioDepXInd = ratioDepXInd, ratioDepYInd = ratioDepYInd,ratioIndDepYIndDepX = ratioIndDepYIndDepX, rep = "rep3")

#########
cond2 = "100nM PLX (3)"
cond1 = "1uM PLX (5)"
filterThreshold <- 0
x = 5; 
y = 3;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x, V2.y=V2.y) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)

cutoff = c(100,200,300,400,500,600,700,800,900,1000)
foldChangecutoffDep = 2.5
foldChangecutoffInd = 1.5

ratioDepXDepY = c()
ratioDepXInd = c()
ratioDepYInd = c()
ratioIndDepYIndDepX = c()

for (i in 1:length(cutoff)) {
  drugDepX <- overlapPLX %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(V2.x>cutoff[i] & V2.y<cutoff[i])
  drugDepY <- overlapPLX %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(V2.y>cutoff[i] & V2.x<cutoff[i])
  drugInd <- overlapPLX %>% filter(foldchange > -foldChangecutoffInd & foldchange < foldChangecutoffInd) %>% filter(V2.x>cutoff[i] & V2.y>cutoff[i])
  
  ratioDepXDepY[i] = nrow(drugDepY)/nrow(drugDepX)
  ratioDepXInd[i] = nrow(drugInd)/nrow(drugDepX)
  ratioDepYInd[i] = nrow(drugInd)/nrow(drugDepY)
  ratioIndDepYIndDepX[i] = (nrow(drugInd) + nrow(drugDepY))/(nrow(drugInd) + nrow(drugDepX))
}
tableReplicate4 = tibble(cutoff = cutoff, ratioDepXDepY = ratioDepXDepY, ratioDepXInd = ratioDepXInd, ratioDepYInd = ratioDepYInd,ratioIndDepYIndDepX = ratioIndDepYIndDepX, rep = "rep4")

tableAll = bind_rows(tableReplicate1, tableReplicate2, tableReplicate3, tableReplicate4)

tableAllMeanSEM = tableAll %>% group_by(cutoff) %>% summarise(meanDepXDepY = mean(ratioDepXDepY), semDepXDepY = sqrt(var(ratioDepXDepY)/length(ratioDepXDepY)), 
                                                              meanIndDepYIndDepX = mean(ratioIndDepYIndDepX), semIndDepYIndDepX = sqrt(var(ratioIndDepYIndDepX)/length(ratioIndDepYIndDepX)))

plot <- ggplot() +
  geom_jitter(aes(x= tableAll$cutoff,y= tableAll$ratioIndDepYIndDepX),alpha=0.6, size = 2.5, shape = 16, width = 15) +
  geom_pointrange(aes(x= tableAllMeanSEM$cutoff, y = tableAllMeanSEM$meanIndDepYIndDepX, ymin=tableAllMeanSEM$meanIndDepYIndDepX-tableAllMeanSEM$semIndDepYIndDepX, ymax=tableAllMeanSEM$meanIndDepYIndDepX+tableAllMeanSEM$semIndDepYIndDepX),size = 0.8, shape = 16)+
  ylim(0, 1.1*max(tableAll$ratioIndDepYIndDepX)) +
  scale_x_continuous(breaks=seq(100,1000,100)) + 
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot = plot, file = paste0(plot3Directory, 'BC03_foldChangeAllratioIndDepYIndDepX','.svg'), width = 10, height = 4)

tableAll400 = tableAll %>% filter(cutoff == 400)
tableAllMeanSEM400 = tableAllMeanSEM %>% filter(cutoff == 400)

plot <- ggplot() +
  geom_jitter(aes(x= tableAll400$cutoff,y= tableAll400$ratioIndDepYIndDepX),alpha=0.6, size = 2.5, shape = 16, width = 12) +
  geom_pointrange(aes(x= tableAllMeanSEM400$cutoff, y = tableAllMeanSEM400$meanIndDepYIndDepX, ymin=tableAllMeanSEM400$meanIndDepYIndDepX-tableAllMeanSEM400$semIndDepYIndDepX, ymax=tableAllMeanSEM400$meanIndDepYIndDepX+tableAllMeanSEM400$semIndDepYIndDepX),size = 0.8, shape = 16)+
  ylim(0, 1.1*max(tableAll400$ratioIndDepYIndDepX)) +
  xlim(370, 430) +
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot = plot, file = paste0(plot3Directory, 'BC03_foldChange400ratioIndDepYIndDepX','.svg'), width = 1.5, height = 4)

