#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/gDNA_BCSeq/20200923gDNA_BCSeq/20200927_BC02/scripts/2929.10.02_BC02_barcodeAnalysis.R
##############################################

library(tidyverse)
library(reshape2)
library(ggridges)

sampleDir <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/gDNA_BarcodeData/BC02/analyzed/'
sampleFolders <- list.files(path=sampleDir,pattern='*d8.txt',recursive=TRUE)
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

#[1] "BC02_A_100nM_PLX_1/starcode/BC02_A_100nM_PLX_1_clusteredBarcodeUMIs_d8.txt"
#[2] "BC02_A_100nM_PLX_2/starcode/BC02_A_100nM_PLX_2_clusteredBarcodeUMIs_d8.txt"
#[3] "BC02_A_1uM_PLX/starcode/BC02_A_1uM_PLX_clusteredBarcodeUMIs_d8.txt"  
#[4] "BC02_B_100nM_PLX/starcode/BC02_B_100nM_PLX_clusteredBarcodeUMIs_d8.txt"    
#[5] "BC02_B_1uM_PLX/starcode/BC02_B_1uM_PLX_clusteredBarcodeUMIs_d8.txt" 

sampleTables <- list(rep("NA",length(sampleFolders)))
for (i in 1:length(sampleFolders)) {
  sampleTable = read.table(paste0(sampleDir, sampleFolders[i]), header = FALSE, sep = "\t")
  sampleTables[[i]] <- sampleTable
}

####Normalizing using spike-in barcodes
spikeBC1 = "TCCAGGTCCTCCTACTTGTACAACACCTTGTACAGCTGCTAGTGGTAGAAGAGGTACAACAACAACACGAGCATCATGAGGATCTACAGCATCAAGAACA"
spikeBC2 = "ACGTTGTGCATGACCTTGATCACCAGCTCGATGTCGAACATCACGAGCTCGTTCTGCATCTGCAAGAACACCTCGTCCTTGAACTGCTCGACGTCCATGA"  
nBC1tot = c(100000,100000,50000,100000,50000) #1-5
nBC2tot = c(20000,20000,10000,20000,10000)    #1-5
nFraction = c(0.32,0.32,0.34,0.32,0.35)    #1-5
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

# ggplot(standardTableAll[c(9:10),], aes(V2, nBC)) + geom_point() +
#   geom_smooth(method = lm, formula=y~0+x, se = FALSE) + xlim(0,4000) + ylim(0,35000)
#   geom_abline(slope=5,intercept=0) 
#   
# ggplot(standardTableAll, aes(c(3002,727), c(32000,6400))) + geom_point()

#CHECK INDEX CONTROL
x = 1; 
y = 2;


indexControlNorm <- full_join(sampleTables[[x]],sampleTables[[y]],by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x*lmr[[x]]) %>% mutate(V2.y=V2.y*lmr[[y]])


#CHECK INDEX CONTROL
indexControlNorm <- full_join(sampleTables[[1]],sampleTables[[2]],by="V1") %>% replace(is.na(.),0) %>% 
  mutate(V2.x=V2.x/sum(V2.x)*10^6) %>% mutate(V2.y=V2.y/sum(V2.y)*10^6)
ggplot() +
  geom_point(aes(x=indexControlNorm$V2.x,y=indexControlNorm$V2.y),alpha=0.1) +
  geom_abline(slope=1,intercept=0) +
  xlab("1uM PLX") + ylab("1uM PLX") +
  theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))





###############################################################################################################################
###########################################DRUG DEPENDENT vs INDEPENDENT LINEAGES##############################################
cond2 = "1uM PLX (3)"
cond1 = "1uM PLX (5)"
filterThreshold <- 0
x = 5; 
y = 3;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x+1, V2.y=V2.y+1) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)

cutoff = 500

drugDepX <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.x>cutoff & V2.y<cutoff)
drugDepY <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>cutoff & V2.x<cutoff)
drugInd <- overlapPLX %>% filter(foldchange > -1.5 & foldchange < 1.5) %>% filter(V2.x>cutoff & V2.y>cutoff)


######For Figure in Paper (Figure 3B)

plot <- ggplot() +
  geom_point(aes(x=overlapPLX$V2.x,y=overlapPLX$V2.y),alpha=0.4, size = 4, shape = 16) +
  geom_point(aes(x=drugInd$V2.x,y=drugInd$V2.y),col="orange", size = 4, shape = 16) +
  geom_point(aes(x=drugDepY$V2.x,y=drugDepY$V2.y),col="springgreen3", size = 4, shape = 16) +
  geom_point(aes(x=drugDepX$V2.x,y=drugDepX$V2.y),col="turquoise3", size = 4, shape = 16) +
  geom_abline(slope=1,intercept=0) +
  xlab(paste0(cond1)) + ylab(paste0(cond2)) +
 # annotate(geom = "text", x = max(overlapPLX$V2.y)/2, y = max(overlapPLX$V2.y)/2, label = expression(italic("  y=x")), hjust = "left") + 
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 2), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot = plot, file = paste0(plotDirectory, 'BC02_independence',cond1,'Vs', cond2,'noAxisNorm_FigureV1500.svg'), width = 8, height = 7)

#####percentages for Venn plot:
percentDepX = nrow(drugDepX)/(nrow(drugDepX) + nrow(drugDepY) + nrow(drugInd))
percentDepY = nrow(drugDepY)/(nrow(drugDepX) + nrow(drugDepY) + nrow(drugInd))
percentdrugInd = nrow(drugInd)/(nrow(drugDepX) + nrow(drugDepY) + nrow(drugInd))

#########################################################################################
##############Analysis to see where the green cells go###################################
#########################################################################################
cond2 = "100nM PLX (1)"
cond1 = "1uM PLX (3)"
filterThreshold <- 0
x = 3; 
y = 1;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x+1, V2.y=V2.y+1) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
drugDepY_1 <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>1500 & V2.x<600) %>% select(V1) %>% unique()

cond2 = "100nM PLX (1)"
cond1 = "1uM PLX (5)"
filterThreshold <- 0
x = 5; 
y = 1;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x+1, V2.y=V2.y+1) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
drugDepY_2 <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>1500 & V2.x<600) %>% select(V1) %>% unique()

cond2 = "100nM PLX (4)"
cond1 = "1uM PLX (3)"
filterThreshold <- 0
x = 3; 
y = 4;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x+1, V2.y=V2.y+1) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
drugDepY_3 <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>1500 & V2.x<600) %>% select(V1) %>% unique()

cond2 = "100nM PLX (4)"
cond1 = "1uM PLX (5)"
filterThreshold <- 0
x = 5; 
y = 4;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x+1, V2.y=V2.y+1) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
drugDepY_4 <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>1500 & V2.x<600) %>% select(V1) %>% unique()

# drugDepY_12 = inner_join(drugDepY_1,drugDepY_2)
# drugDepY_34 = inner_join(drugDepY_3,drugDepY_4)

drugDepY_All <- bind_rows(drugDepY_1,drugDepY_2,drugDepY_3,drugDepY_4) %>% select(V1) %>% unique()
# drugDepY_AllB <- bind_rows(drugDepY_12,drugDepY_34) %>% unique()

cond2 = "100nM PLX (1) only-Dependent"
cond1 = "100nM PLX (4) only-Dependent"
filterThreshold <- 0
x = 4; 
y = 1;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x+1, V2.y=V2.y+1) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)

drugDepX <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.x>1500 & V2.y<600)
drugDepY <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>1500 & V2.x<600)
drugInd <- overlapPLX %>% filter(foldchange > -1.5 & foldchange < 1.5) %>% filter(V2.x>400 & V2.y>400)
drugDepY_All_joined = inner_join(overlapPLX,drugDepY_All, by = "V1")
# drugDepY_All_joinedB = inner_join(overlapPLX,drugDepY_AllB, by = "V1")

drugDepX1 <- drugDepY_All_joined %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.x>1500 & V2.y<600)
drugDepY1 <- drugDepY_All_joined %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>1500 & V2.x<600)
drugInd1 <- drugDepY_All_joined %>% filter(foldchange > -1.5 & foldchange < 1.5) %>% filter(V2.x>400 & V2.y>400)

plot <- ggplot() +
  geom_point(aes(x=overlapPLX$V2.x,y=overlapPLX$V2.y),alpha=0.2) +
  geom_point(aes(x=drugDepY_All_joined$V2.x,y=drugDepY_All_joined$V2.y),col="orange") +
  geom_point(aes(x=drugDepY1$V2.x,y=drugDepY1$V2.y),col="springgreen3") +
  geom_point(aes(x=drugDepX1$V2.x,y=drugDepX1$V2.y),col="turquoise3") +
  geom_abline(slope=1,intercept=0) +
  xlab(paste0(cond1)) + ylab(paste0(cond2))+
  annotate(geom = "text", x = max(overlapPLX$V2.y)/2, y = max(overlapPLX$V2.y)/2, label = expression(italic("  y=x")), hjust = "left") + 
  theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))
#ggsave(plot = plot, file = paste0(plotDirectory, 'independence',cond1,'Vs', cond2,'.pdf'), width = 8, height = 7)


plot <-  ggplot() +
  geom_point(aes(x=drugDepY_All_joined$V2.x,y=drugDepY_All_joined$V2.y),alpha=0.2) +
  geom_point(aes(x=drugInd$V2.x,y=drugInd$V2.y),col="orange") +
  geom_point(aes(x=drugDepY$V2.x,y=drugDepY$V2.y),col="springgreen3") +
  geom_point(aes(x=drugDepX$V2.x,y=drugDepX$V2.y),col="turquoise3") +
  geom_abline(slope=1,intercept=0) +
  xlab(paste0(cond1)) + ylab(paste0(cond2))+ 
  annotate(geom = "text", x = max(overlapPLX$V2.y)/2, y = max(overlapPLX$V2.y)/2, label = expression(italic("  y=x")), hjust = "left") + 
  theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))
#ggsave(plot = plot, file = paste0(plotDirectory, 'independence',cond1,'Vs', cond2,'.pdf'), width = 8, height = 7)

#########################################################################################
##############Unique number of barcodes per conditions###################################
#########################################################################################
#[1] "BC02_A_100nM_PLX_1/starcode/BC02_A_100nM_PLX_1_clusteredBarcodeUMIs_d8.txt"
#[2] "BC02_A_100nM_PLX_2/starcode/BC02_A_100nM_PLX_2_clusteredBarcodeUMIs_d8.txt"
#[3] "BC02_A_1uM_PLX/starcode/BC02_A_1uM_PLX_clusteredBarcodeUMIs_d8.txt"        
#[4] "BC02_B_100nM_PLX/starcode/BC02_B_100nM_PLX_clusteredBarcodeUMIs_d8.txt"    
#[5] "BC02_B_1uM_PLX/starcode/BC02_B_1uM_PLX_clusteredBarcodeUMIs_d8.txt" 

nFraction = c(0.32,0.32,0.34,0.32,0.35)    #1-5

BC02_100nMPLX_A1 = length(unique(sampleTables[[1]]$V1))/nFraction[1]
BC02_100nMPLX_A2 = length(unique(sampleTables[[2]]$V1))/nFraction[2]
BC02_1uMPLX_A1 = length(unique(sampleTables[[3]]$V1))/nFraction[3]
BC02_100nMPLX_B1 = length(unique(sampleTables[[4]]$V1))/nFraction[4]
BC02_1uMPLX_B1 = length(unique(sampleTables[[5]]$V1))/nFraction[5]



#########################################################################################
###########################################DRUG DEPENDENT vs INDEPENDENT LINEAGES##############################################
###########################################OverCutoff##############################################

# (2,3) 100nM PLX vs 1uM PLX
cond2 = "100nM PLX (1)"
cond1 = "1uM PLX (3)"
filterThreshold <- 0
x = 3; 
y = 1;
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
}
tableReplicate1 = tibble(cutoff = cutoff, ratioDepXDepY = ratioDepXDepY, ratioDepXInd = ratioDepXInd, ratioDepYInd = ratioDepYInd,ratioIndDepYIndDepX = ratioIndDepYIndDepX, rep = "rep1")

# (2,3) 100nM PLX vs 1uM PLX
cond2 = "100nM PLX (4)"
cond1 = "1uM PLX (3)"
filterThreshold <- 0
x = 3; 
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
tableReplicate3 = tibble(cutoff = cutoff, ratioDepXDepY = ratioDepXDepY, ratioDepXInd = ratioDepXInd, ratioDepYInd = ratioDepYInd,ratioIndDepYIndDepX = ratioIndDepYIndDepX, rep = "rep3")

#########
cond2 = "100nM PLX (1)"
cond1 = "1uM PLX (5)"
filterThreshold <- 0
x = 5; 
y = 1;
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
  ylim(0, 1.2*max(tableAll$ratioIndDepYIndDepX)) +
  scale_x_continuous(breaks=seq(100,1000,100)) + 
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot = plot, file = paste0(plot3Directory, 'BC02_foldChangeAllratioIndDepYIndDepX','.svg'), width = 10, height = 4)


########OLD
# ###########################################################################################################
# #################### PAIRWISE HERITABILITY ################################################################
# # HERITABILITY- (1,3) 100nM PLX vs 1uM PLX
# cond1 = "100nM PLX (1)"
# cond2 = "1uM PLX (3)"
# filterThreshold <- 0
# x = 1; 
# y = 3;
# a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
# b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
# herit1uMPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0)
# scatterPlot<- ggplot() +
#   geom_point(aes(x=herit1uMPLX$V2.x,y=herit1uMPLX$V2.y),alpha=0.3) +
#   geom_abline(slope=1,intercept=0) +
#   xlab("100nM PLX") + ylab("1uM PLX") +
#   annotate(geom = "text", x = max(herit1uMPLX$V2.x)/2, y = max(herit1uMPLX$V2.x)/2, label = expression(italic(" y=x")), hjust = "left") + 
#   theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))
# ggsave(plot = scatterPlot, file = paste0(plotDirectory, 'scatter',cond1,'Vs', cond2,'.pdf'), width = 8, height = 7)
# 
# overlapDMSOTable <- makeOverlapTable(a,b) + xlab(paste0("Top n Barcodes in ",cond1)) +
#   ylab(paste0("Top n Barcodes in ",cond2))
# 
# ggsave(plot = overlapDMSOTable, file = paste0(plotDirectory, 'heritability',cond1,'Vs', cond2,'.pdf'), width = 7, height = 7)
# 
# # HERITABILITY- (2,3) 100nM PLX vs 1uM PLX
# cond1 = "100nM PLX (2)"
# cond2 = "1uM PLX (3)"
# x = 2; 
# y = 3;
# a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
# b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
# herit1uMPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0)
# scatterPlot<- ggplot() +
#   geom_point(aes(x=herit1uMPLX$V2.x,y=herit1uMPLX$V2.y),alpha=0.3) +
#   geom_abline(slope=1,intercept=0) +
#   xlab("100nM PLX") + ylab("1uM PLX") +
#   annotate(geom = "text", x = max(herit1uMPLX$V2.x)/2, y = max(herit1uMPLX$V2.x)/2, label = expression(italic(" y=x")), hjust = "left") + 
#   theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))
# ggsave(plot = scatterPlot, file = paste0(plotDirectory, 'scatter',cond1,'Vs', cond2,'.pdf'), width = 8, height = 7)
# 
# overlapDMSOTable <- makeOverlapTable(a,b) + xlab(paste0("Top n Barcodes in ",cond1)) +
#   ylab(paste0("Top n Barcodes in ",cond2))
# ggsave(plot = overlapDMSOTable, file = paste0(plotDirectory, 'heritability',cond1,'Vs', cond2,'.pdf'), width = 7, height = 7)
# 
# #CHECK HERITABILITY- (2,4) 100nM PLX vs 100nM PLX
# cond1 = "100nM PLX (2)"
# cond2 = "100nM PLX (4)"
# x = 2; 
# y = 4;
# a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
# b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
# herit1uMPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0)
# scatterPlot<- ggplot() +
#   geom_point(aes(x=herit1uMPLX$V2.x,y=herit1uMPLX$V2.y),alpha=0.3) +
#   geom_abline(slope=1,intercept=0) +
#   xlab("100nM PLX") + ylab("1uM PLX") +
#   annotate(geom = "text", x = max(herit1uMPLX$V2.x)/2, y = max(herit1uMPLX$V2.x)/2, label = expression(italic(" y=x")), hjust = "left") + 
#   theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))
# ggsave(plot = scatterPlot, file = paste0(plotDirectory, 'scatter',cond1,'Vs', cond2,'.pdf'), width = 8, height = 7)
# 
# overlapDMSOTable <- makeOverlapTable(a,b) + xlab(paste0("Top n Barcodes in ",cond1)) +
#   ylab(paste0("Top n Barcodes in ",cond2))
# ggsave(plot = overlapDMSOTable, file = paste0(plotDirectory, 'heritability',cond1,'Vs', cond2,'.pdf'), width = 7, height = 7)
# 
# #CHECK HERITABILITY- (2,5) 100nM PLX vs 1uM PLX
# cond1 = "100nM PLX (2)"
# cond2 = "1uM PLX (5)"
# x = 2; 
# y = 5;
# a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
# b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
# herit1uMPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0)
# scatterPlot<- ggplot() +
#   geom_point(aes(x=herit1uMPLX$V2.x,y=herit1uMPLX$V2.y),alpha=0.3) +
#   geom_abline(slope=1,intercept=0) +
#   xlab("100nM PLX") + ylab("1uM PLX") +
#   annotate(geom = "text", x = max(herit1uMPLX$V2.x)/2, y = max(herit1uMPLX$V2.x)/2, label = expression(italic(" y=x")), hjust = "left") + 
#   theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))
# ggsave(plot = scatterPlot, file = paste0(plotDirectory, 'scatter',cond1,'Vs', cond2,'.pdf'), width = 8, height = 7)
# 
# overlapDMSOTable <- makeOverlapTable(a,b) + xlab(paste0("Top n Barcodes in ",cond1)) +
#   ylab(paste0("Top n Barcodes in ",cond2))
# ggsave(plot = overlapDMSOTable, file = paste0(plotDirectory, 'heritability',cond1,'Vs', cond2,'.pdf'), width = 7, height = 7)
# 
# #CHECK HERITABILITY- (3,5) 1uM PLX vs 1uM PLX
# cond1 = "1uM PLX (3)"
# cond2 = "1uM PLX (5)"
# x = 3; 
# y = 5;
# a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
# b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
# herit1uMPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0)
# scatterPlot<- ggplot() +
#   geom_point(aes(x=herit1uMPLX$V2.x,y=herit1uMPLX$V2.y),alpha=0.3) +
#   geom_abline(slope=1,intercept=0) +
#   xlab("100nM PLX") + ylab("1uM PLX") +
#   annotate(geom = "text", x = max(herit1uMPLX$V2.x)/2, y = max(herit1uMPLX$V2.x)/2, label = expression(italic(" y=x")), hjust = "left") + 
#   theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))
# ggsave(plot = scatterPlot, file = paste0(plotDirectory, 'scatter',cond1,'Vs', cond2,'.pdf'), width = 8, height = 7)
# 
# overlapDMSOTable <- makeOverlapTable(a,b) + xlab(paste0("Top n Barcodes in ",cond1)) +
#   ylab(paste0("Top n Barcodes in ",cond2))
# ggsave(plot = overlapDMSOTable, file = paste0(plotDirectory, 'heritability',cond1,'Vs', cond2,'.pdf'), width = 7, height = 7)
# 
# # HERITABILITY- (1,2) 100nM PLX vs 1uM PLX
# cond1 = "1uM PLX (1)"
# cond2 = "1uM PLX (2)"
# filterThreshold <- 0
# x = 1; 
# y = 2;
# a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
# b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
# herit1uMPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0)
# scatterPlot<- ggplot() +
#   geom_point(aes(x=herit1uMPLX$V2.x,y=herit1uMPLX$V2.y),alpha=0.3) +
#   geom_abline(slope=1,intercept=0) +
#   xlab("100nM PLX") + ylab("1uM PLX") +
#   annotate(geom = "text", x = max(herit1uMPLX$V2.x)/2, y = max(herit1uMPLX$V2.x)/2, label = expression(italic(" y=x")), hjust = "left") + 
#   theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))
# ggsave(plot = scatterPlot, file = paste0(plotDirectory, 'scatter',cond1,'Vs', cond2,'.pdf'), width = 8, height = 7)
# 
# overlapDMSOTable <- makeOverlapTable(a,b) + xlab(paste0("Top n Barcodes in ",cond1)) +
#   ylab(paste0("Top n Barcodes in ",cond2))
# 
# ggsave(plot = overlapDMSOTable, file = paste0(plotDirectory, 'heritability',cond1,'Vs', cond2,'.pdf'), width = 7, height = 7)
# 
# ###############################################################################################################################
# ###########################################DRUG DEPENDENT vs INDEPENDENT LINEAGES##############################################
# # (1,3) 100nM PLX vs 1uM PLX
# cond2 = "100nM PLX (1)"
# cond1 = "1uM PLX (3)"
# filterThreshold <- 0
# x = 3; 
# y = 1;
# a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
# b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
# overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x, V2.y=V2.y) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
# 
# drugDepX <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.x>1500 & V2.y<600)
# drugDepY <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>1500 & V2.x<600)
# drugInd <- overlapPLX %>% filter(foldchange > -1.5 & foldchange < 1.5) %>% filter(V2.x>400 & V2.y>400)
# 
# plot <- ggplot() +
#   geom_point(aes(x=overlapPLX$V2.x,y=overlapPLX$V2.y),alpha=0.2) +
#   geom_point(aes(x=drugInd$V2.x,y=drugInd$V2.y),col="orange") +
#   geom_point(aes(x=drugDepY$V2.x,y=drugDepY$V2.y),col="springgreen3") +
#   geom_point(aes(x=drugDepX$V2.x,y=drugDepX$V2.y),col="turquoise3") +
#   geom_abline(slope=1,intercept=0) +
#   xlab(paste0(cond1)) + ylab(paste0(cond2)) +
#   annotate(geom = "text", x = max(overlapPLX$V2.y)/2, y = max(overlapPLX$V2.y)/2, label = expression(italic("  y=x")), hjust = "left") + 
#   theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))
# ggsave(plot = plot, file = paste0(plotDirectory, 'independence',cond1,'Vs', cond2,'noAxisNorm.pdf'), width = 8, height = 7)
# 
# ######
# # (1,3) 100nM PLX vs 1uM PLX
# cond2 = "100nM PLX (4)"
# cond1 = "1uM PLX (5)"
# filterThreshold <- 0
# x = 5; 
# y = 4;
# a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
# b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
# overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x+1, V2.y=V2.y+1) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)
# 
# drugDepX <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.x>1500 & V2.y<600)
# drugDepY <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>1500 & V2.x<600)
# drugInd <- overlapPLX %>% filter(foldchange > -1.5 & foldchange < 1.5) %>% filter(V2.x>400 & V2.y>400)
# 
# plot <- ggplot() +
#   geom_point(aes(x=overlapPLX$V2.x,y=overlapPLX$V2.y),alpha=0.2) +
#   geom_point(aes(x=drugInd$V2.x,y=drugInd$V2.y),col="orange") +
#   geom_point(aes(x=drugDepY$V2.x,y=drugDepY$V2.y),col="springgreen3") +
#   geom_point(aes(x=drugDepX$V2.x,y=drugDepX$V2.y),col="turquoise3") +
#   geom_abline(slope=1,intercept=0) +
#   xlab(paste0(cond1)) + ylab(paste0(cond2)) +
#   annotate(geom = "text", x = max(overlapPLX$V2.y)/2, y = max(overlapPLX$V2.y)/2, label = expression(italic("  y=x")), hjust = "left") + 
#   theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))
# ggsave(plot = plot, file = paste0(plotDirectory, 'independence',cond1,'Vs', cond2,'noAxisNorm.pdf'), width = 8, height = 7)
# 
