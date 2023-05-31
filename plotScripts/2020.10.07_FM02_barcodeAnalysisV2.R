#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/gDNA_BCSeq/20200923gDNA_BCSeq/20201112_FM02/scripts/2020.10.07_FM02_barcodeAnalysis.R
##############################################

library(tidyverse)
library(reshape2)
library(ggridges)
library(spgs)

sampleDir <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/gDNA_BarcodeData/FM02/analyzed/'
sampleFolders <- list.files(path=sampleDir,pattern='*d8.txt',recursive=TRUE)
plotDirectory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/plotData/analysisPlots/'
dataDirectory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/gDNA_BarcodeData/FM02/dataFor10X/'
data2Directory <-"/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM02/"

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

#[1] "FM02_100nM_PLX_5nM_T/starcode/FM02_100nM_PLX_5nM_T_clusteredBarcodeUMIs_d8.txt" 
#[2] "FM02_100nM_PLX/starcode/FM02_100nM_PLX_clusteredBarcodeUMIs_d8.txt"            
#[3] "FM02_1uM_PLX/starcode/FM02_1uM_PLX_clusteredBarcodeUMIs_d8.txt"                 
#[4] "FM02_5nM_T/starcode/FM02_5nM_T_clusteredBarcodeUMIs_d8.txt"    

sampleTables <- list(rep("NA",length(sampleFolders)))
for (i in 1:length(sampleFolders)) {
  sampleTable = read.table(paste0(sampleDir, sampleFolders[i]), header = FALSE, sep = "\t")
  sampleTables[[i]] <- sampleTable
}

####Normalizing using spike-in barcodes
spikeBC1 = "TCCAGGTCCTCCTACTTGTACAACACCTTGTACAGCTGCTAGTGGTAGAAGAGGTACAACAACAACACGAGCATCATGAGGATCTACAGCATCAAGAACA"
spikeBC2 = "ACGTTGTGCATGACCTTGATCACCAGCTCGATGTCGAACATCACGAGCTCGTTCTGCATCTGCAAGAACACCTCGTCCTTGAACTGCTCGACGTCCATGA"  
nBC1tot = c(10000,100000,10000,10000) #1-4
nBC2tot = c(2000,20000,2000,2000)    #1-4
nFraction = c(0.76,0.43,0.48,0.46)    #1-4
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
# (2,3) 100nM PLX vs 1uM PLX
cond2 = "100nM PLX (2)"
cond1 = "1uM PLX (3)"
filterThreshold <- 0
x = 3; 
y = 2;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[x]])
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2*lmr[[y]])
overlapPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0) %>% mutate(V2.x=V2.x, V2.y=V2.y) %>% mutate(foldchange = log2((V2.y+1)/(V2.x+1))) %>% mutate(combo = V2.x + V2.y)

cutoff = 400

drugDepX <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.x>cutoff & V2.y<cutoff)
drugDepY <- overlapPLX %>% filter(foldchange < -2.5 | foldchange > 2.5) %>% filter(V2.y>cutoff & V2.x<cutoff)
drugInd <- overlapPLX %>% filter(foldchange > -1.5 & foldchange < 1.5) %>% filter(V2.x>cutoff & V2.y>cutoff)


plot <- ggplot() +
  geom_point(aes(x=overlapPLX$V2.x,y=overlapPLX$V2.y),alpha=0.4, size = 4, shape = 16) +
  geom_point(aes(x=drugInd$V2.x,y=drugInd$V2.y),col="orange", size = 4, shape = 16) +
  geom_point(aes(x=drugDepY$V2.x,y=drugDepY$V2.y),col="springgreen3", size = 4, shape = 16) +
  geom_point(aes(x=drugDepX$V2.x,y=drugDepX$V2.y),col="turquoise3", size = 4, shape = 16) +
  geom_abline(slope=1,intercept=0) +
  geom_abline(slope=1,intercept=0) +
  xlab(paste0(cond1)) + ylab(paste0(cond2)) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 2), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(plot = plot, file = paste0(plotDirectory, 'FM02_independence',cond1,'Vs', cond2,'noAxisNorm.svg'), width = 8, height = 7)


########################################################################################################################################
###########################################GENERATING BARCODES FOR 10X BARCODE COMPARISONS##############################################
########################################################################################################################################
lowDoseDependentBarcodes <- drugDepY$V1
lowDoseDependentBarcodesReverseComplement = as_tibble(reverseComplement(lowDoseDependentBarcodes, content="dna",case="as is"))
lowDoseDependentBarcodes <- as_tibble(drugDepY$V1)
write.table(lowDoseDependentBarcodesReverseComplement, file=paste0(dataDirectory,'lowDoseDependentBarcodesReverseComplement.csv'), sep=",", row.names = FALSE, col.names = FALSE)
write.table(lowDoseDependentBarcodes, file=paste0(dataDirectory,'lowDoseDependentBarcodes.csv'), sep=",", row.names = FALSE, col.names = FALSE)


########################################################################################################################################
###########################################GENERATING BARCODES FOR 10X BARCODE COMPARISONS ___PRESENT IN BOTH##############################################
########################################################################################################################################
bothLowHigh <- drugInd$V1
bothLowHighReverseComplement = as_tibble(reverseComplement(bothLowHigh, content="dna",case="as is"))
bothLowHigh <- as_tibble(drugInd$V1)
write.table(bothLowHighReverseComplement, file=paste0(dataDirectory,'bothLowHighReverseComplement.csv'), sep=",", row.names = FALSE, col.names = FALSE)
write.table(bothLowHigh, file=paste0(dataDirectory,'bothLowHigh.csv'), sep=",", row.names = FALSE, col.names = FALSE)




