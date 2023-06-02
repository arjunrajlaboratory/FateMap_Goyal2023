#####################NOTE#####################
#Cleaned up version of V2: /Users/yogesh/Dropbox (RajLab)/FateMap/data/gDNA_BCSeq/20210211_Transduction_BC04/20210211_Transductions/scripts/2021.02.12_Transduction_barcodeAnalysis.R
##############################################

library(tidyverse)
library(reshape2)
library(ggridges)


sampleDir <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/gDNA_BarcodeData/BC04_Transductions/analyzed/'
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

#[1] "transduction_1_2/starcode/transduction_1_2_clusteredBarcodeUMIs_d8.txt" 
#[2] "transduction_3_5/starcode/transduction_3_5_clusteredBarcodeUMIs_d8.txt"            
#[3] "transduction_8_9/starcode/transduction_8_9_clusteredBarcodeUMIs_d8.txt"                 


sampleTables <- list(rep("NA",length(sampleFolders)))
for (i in 1:length(sampleFolders)) {
  sampleTable = read.table(paste0(sampleDir, sampleFolders[i]), header = FALSE, sep = "\t")
  sampleTables[[i]] <- sampleTable
}


###########################################################################################################
#################### PAIRWISE HERITABILITY ################################################################
# HERITABILITY- Transductions_1_2 vs Transductions_3_5
cond1 = "Transductions_1_2"
cond2 = "Transductions_3_5"
filterThreshold <- 1
x = 1;
y = 2;
z = 3;
a <- sampleTables[[x]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*(10^6))
b <- sampleTables[[y]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*(10^6))
c <- sampleTables[[z]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*(10^6))

herit1uMPLX <- full_join(a,b,by="V1") %>% replace(is.na(.),0)
test <- inner_join(a,b,by="V1") %>% replace(is.na(.),0)

estimate1 = nrow(a)/(nrow(inner_join(a,b,by="V1"))/(nrow(b)))
estimate2 = nrow(a)/(nrow(inner_join(a,c,by="V1"))/(nrow(c)))
estimate3 = nrow(b)/(nrow(inner_join(b,c,by="V1"))/(nrow(c)))

a1 = sample_n(a,4500)
b1 = sample_n(b,7500)
herit1uMPLX <- full_join(a1,b1,by="V1") %>% replace(is.na(.),0)
test <- inner_join(a1,b1,by="V1") %>% replace(is.na(.),0)

scatterPlot = ggplot() +
  geom_point(aes(x=herit1uMPLX$V2.x,y=herit1uMPLX$V2.y),alpha=0.3) +
  geom_abline(slope=1,intercept=0) +
  xlab(paste0(cond1)) + ylab(paste0(cond2)) +
  annotate(geom = "text", x = max(herit1uMPLX$V2.x)/2, y = max(herit1uMPLX$V2.x)/2, label = expression(italic(" y=x")), hjust = "left") + 
  theme_classic((base_size = 18)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 1), text=element_text(family="Helvetica"))
ggsave(plot = scatterPlot, file = paste0(plotDirectory, 'scatter',cond1,'Vs', cond2,'.pdf'), width = 8, height = 7)