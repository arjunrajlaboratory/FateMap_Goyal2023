#####################NOTE#####################
#Cleaned up version of V1: /Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2019_FM01/SubmissionScripts/seurat_analysis/memory_s1s2_integration_GlobalAnalysisV1.R
##############################################

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RANN)
library(tripack)
library(reshape2)
library(svglite)

home1Directory <-("/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_SeuratObject/FM01/")
home2Directory <- '/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/10X_BarcodeData/FM01/'
plot5Directory <- "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/plotData/analysisPlots/"

pcaCoordinates = as_tibble(read.table(file = paste0(home1Directory, "pcaCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))

####functions

create_lpr_theme <- function(){
  lpr_theme <- ggplot2::theme_bw() + ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 32),
                   plot.title = ggplot2::element_text(face = "bold",
                                                      hjust = 0.5, size = ggplot2::rel(0.5)), axis.text = ggplot2::element_text(size = ggplot2::rel(0.8),
                                                                                                                                color = "black", face = "bold"), axis.title = ggplot2::element_text(size = ggplot2::rel(0.6),
                                                                                                                                                                                                    face = "bold"), legend.title = ggplot2::element_blank(),
                   legend.position = "bottom", axis.line = element_line(size = 2),
                   axis.ticks = element_line(size = 2), strip.background = element_rect(size = 2),
                   strip.text = element_text(size = ggplot2::rel(0.7),
                                             face = "bold"))
  return(lpr_theme)
}

###### Get lower triangle of the correlation matrix
get_lower_tri<-function(neightborMatrix2){
  neightborMatrix2[upper.tri(neightborMatrix2)] <- NA
  return(neightborMatrix2)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(neightborMatrix){
  neightborMatrix[lower.tri(neightborMatrix)]<- NA
  return(neightborMatrix)
}

#######Barcode entries for Sample 1 and Sample 2
#################################################
barcode50 = as_tibble(read.table(paste0(home2Directory, 'stepFourStarcodeShavedReads50.txt'), stringsAsFactors=F, header = T))
umiCut = 4 # Minimum UMI cutoff for reliable analysis. Can be played around with between 3 to 5.
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.  

upToLineageCounting = barcode50 %>% 
  filter(sampleNum %in% c("1","2")) %>%
  group_by(cellID, BC50StarcodeD8, sampleNum) %>%
  summarise(nUMI = length(sampleNum)) %>%
  filter(nUMI >= umiCut) %>% 
  group_by(cellID) %>%
  mutate(nLineages = length(cellID)) 

#making the nomenclature consistent with the Matrix from data integration
upToLineageCounting$sampleNum = sub("^", "S", upToLineageCounting$sampleNum) 
#####
###taking only single cellid-->barcode mappings
linCountTooverlaps = upToLineageCounting %>%
  ungroup() %>%
  filter(nLineages <= linCut) %>%
  unique()

##################################################################################################
jointPCA = inner_join(linCountTooverlaps, pcaCoordinates, by = c("cellID", "sampleNum")) %>% select(-nLineages)
#BarcodesRename = jointPCA %>% select(BC50StarcodeD8) %>% unique() 
##################################################################################################
jointPCA1Barcodes = jointPCA %>% filter(sampleNum=="S1") %>% select(BC50StarcodeD8) %>% unique()
jointPCA2Barcodes = jointPCA %>% filter(sampleNum=="S2")  %>% select(BC50StarcodeD8) %>% unique()

jointBarcodesOnlyBoth = inner_join(jointPCA2Barcodes,jointPCA1Barcodes, by = "BC50StarcodeD8") 
jointBarcodesOnlyBoth = jointBarcodesOnlyBoth %>% mutate(barcodeName = c(1:nrow(jointBarcodesOnlyBoth)))
jointBarcodesOnlyBoth$barcodeName = sub("^", "B", jointBarcodesOnlyBoth$barcodeName)

jointBarcodesOnlyBothList = jointBarcodesOnlyBoth$BC50StarcodeD8
jointPCAOnlyBoth = jointPCA %>% filter(BC50StarcodeD8 %in% jointBarcodesOnlyBothList)
finalPCAJoint = inner_join(jointPCAOnlyBoth,jointBarcodesOnlyBoth, by = "BC50StarcodeD8") %>% select(-cellID, -BC50StarcodeD8, -nUMI)  ####table of PCs with SampleNum (1), BarcodeName (52)

finalPCAS1 = finalPCAJoint %>% filter(sampleNum=="S1")
finalPCAS2 = finalPCAJoint %>% filter(sampleNum=="S2")

finalPCAS1Big = finalPCAS1 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony >2) %>% select(-nColony)
finalPCAS2Big = finalPCAS2 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony >2) %>% select(-nColony)

commonS1S2Big = inner_join(finalPCAS1Big, finalPCAS2Big)
commonS1S2Big = commonS1S2Big$barcodeName
finalPCAJointBig = finalPCAJoint %>% filter(barcodeName %in% commonS1S2Big)

#####Calculating values for All Siblings ########
meanNeighbor = tibble(withSelf = numeric(),
                      withOther = numeric(),
                      BarcodeName = character());
medianNeighbor = tibble(withSelf = numeric(),
                      withOther = numeric(),
                      BarcodeName = character());
for (i in c(1:length(unique(finalPCAJointBig$barcodeName)))) {
  barcode = unique(finalPCAJointBig$barcodeName)
  BarcodeBx = finalPCAJointBig %>% filter(barcodeName ==barcode[i]) %>% select(-barcodeName)
  BarcodeRef = BarcodeBx %>% mutate(num = c(1:nrow(BarcodeBx)))
  
  S1index = BarcodeRef %>% filter(sampleNum == "S1") %>% select(num)
  S2index = BarcodeRef %>% filter(sampleNum == "S2") %>% select(num)
  ####
  knnPCA = nn2(BarcodeBx[,2:51], BarcodeBx[,2:51])
  neighborKNN = as_tibble(knnPCA[[1]]) ####gets the neighbor
  nS1 = c();
  nS2 = c();
  fraction = (nrow(S2index)/nrow(S1index))
  for (j in c(1:nrow(BarcodeBx))) {
    nS1[j] =  (sum(neighborKNN[j,2:10] %in% S1index$num))*(nrow(S2index)/nrow(S1index))  #this part is changed
    nS2[j] =  sum(neighborKNN[j,2:10] %in% S2index$num)
  }
  neighbor = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(BarcodeBx)))
  neighbor = inner_join(neighbor, BarcodeRef, by = "num")
  neighborS1 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(BarcodeBx))) %>% filter(num %in% S1index$num)
  neighborS2 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(BarcodeBx))) %>% filter(num %in% S2index$num)
  S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(isSelf = "withSelf")  ###if want to include individual   S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(relation = "S1-S1", isSelf = "withSelf")
  S2S2 = as_tibble(neighborS2$fractionS2) %>% mutate(isSelf = "withSelf")
  S1S2 = as_tibble((neighborS1$fractionS2)) %>% mutate(isSelf ="withOther")
  S2S1 = as_tibble((neighborS2$fractionS1)) %>% mutate(isSelf = "withOther")
  toPlot = bind_rows(S1S1,S2S2,S1S2,S2S1) %>% rename(fraction = value) %>% mutate(isSelf = factor(isSelf, levels=c("withSelf","withOther")))
  withSelf = toPlot %>% filter(isSelf == "withSelf") %>% select(-isSelf)
  withOther = toPlot %>% filter(isSelf == "withOther") %>% select(-isSelf)
  meanNeighbor = meanNeighbor %>% add_row(withSelf = mean(withSelf$fraction), 
                                          withOther = mean(withOther$fraction),
                                          BarcodeName = barcode[i])
  medianNeighbor = medianNeighbor %>% add_row(withSelf = median(withSelf$fraction), 
                                          withOther = median(withOther$fraction),
                                          BarcodeName = barcode[i])
}

meanNeighbor = meanNeighbor %>% mutate(ratio = withOther/withSelf)
medianNeighbor = medianNeighbor %>% mutate(ratio = withOther/withSelf)

##############################Figure 2E#####################################
#####Calculating values for All Barcodes from within sibling cohort ########
################################################################################################
################################################################################################
neightborMatrix = matrix(, nrow = length(unique(finalPCAJointBig$barcodeName)), ncol = length(unique(finalPCAJointBig$barcodeName)))
for (i in c(1:length(unique(finalPCAJointBig$barcodeName)))) {
  for (j in c(1:length(unique(finalPCAJointBig$barcodeName)))) {
    barcode = unique(finalPCAJointBig$barcodeName)
    iBarcodeBx = finalPCAJointBig %>% filter(sampleNum == "S1", barcodeName ==barcode[i]) %>% select(-barcodeName)
    jBarcodeBx = finalPCAJointBig %>% filter(sampleNum == "S2", barcodeName ==barcode[j]) %>% select(-barcodeName)
    ijBarcodeBx = bind_rows(iBarcodeBx,jBarcodeBx)
    BarcodeRef = ijBarcodeBx %>% mutate(num = c(1:nrow(ijBarcodeBx)))
    S1index = BarcodeRef %>% filter(sampleNum == "S1") %>% select(num)
    S2index = BarcodeRef %>% filter(sampleNum == "S2") %>% select(num)
    
    knnPCA = nn2(ijBarcodeBx[,2:51], ijBarcodeBx[,2:51])
    neighborKNN = as_tibble(knnPCA[[1]]) ####gets the neighbor
    nS1 = c();
    nS2 = c();
    fraction2 = (nrow(S2index)/nrow(S1index))
    for (k in c(1:nrow(ijBarcodeBx))) {
      nS1[k] =  (sum(neighborKNN[k,2:10] %in% S1index$num))*(nrow(S2index)/nrow(S1index))  #this part is changed
      nS2[k] =  sum(neighborKNN[k,2:10] %in% S2index$num)
    }
    neighbor = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx)))
    neighbor = inner_join(neighbor, BarcodeRef, by = "num")
    neighborS1 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx))) %>% filter(num %in% S1index$num)
    neighborS2 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx))) %>% filter(num %in% S2index$num)
    S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(isSelf = "withSelf")  ###if want to include individual   S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(relation = "S1-S1", isSelf = "withSelf")
    S2S2 = as_tibble(neighborS2$fractionS2) %>% mutate(isSelf = "withSelf")
    S1S2 = as_tibble((neighborS1$fractionS2)) %>% mutate(isSelf ="withOther")
    S2S1 = as_tibble((neighborS2$fractionS1)) %>% mutate(isSelf = "withOther")
    toPlot = bind_rows(S1S1,S2S2,S1S2,S2S1) %>% rename(fraction = value) %>% mutate(isSelf = factor(isSelf, levels=c("withSelf","withOther")))
    withSelf = toPlot %>% filter(isSelf == "withSelf") %>% select(-isSelf)
    withOther = toPlot %>% filter(isSelf == "withOther") %>% select(-isSelf)
    neightborMatrix[i,j]= mean(withOther$fraction)/mean(withSelf$fraction)
  }
}

diagonal = diag(neightborMatrix)
diagonal = diagonal[1:12]
diagonalMatrix = diag(diagonal)
diagonalMatrix[col(diagonalMatrix)!=row(diagonalMatrix)] = NA

neightborMatrixRevised = neightborMatrix[1:12,1:12]

rowNames = c(1:12)
rowNames = sub("^", "Twin ", rowNames)
colNames = rowNames
rownames(diagonalMatrix) = rowNames
colnames(diagonalMatrix) = colNames
melted_Matrix <- melt(diagonalMatrix, na.rm = TRUE)
melted_Matrix$value <- round(melted_Matrix$value,2)

plot = ggplot(data = melted_Matrix, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "black", 
                       limit = c(0,1), 
                       breaks = c(0,0.5,1),
                       labels = c(0,0.5,1),
                       space = "Lab", 
                       name="Mixing\nCoefficient") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
  theme_classic((base_size = 26)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(colour = 'black', size = 1.5),
    axis.text.x = element_text(angle = 90),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(colour = "black", size = 1.5),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.8),
    legend.direction = "horizontal",
    text=element_text(family="Helvetica"))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  scale_y_discrete(position = "right")+
  theme(legend.position = "none")
ggsave(plot, file = paste0(plot5Directory, 'FM01_mixingCoefficientMatrixTwins.svg'), width =8, height = 8)


rowNames = c(1:12)
rowNames = sub("^", "Twin ", rowNames)
colNames = rowNames
rownames(neightborMatrixRevised) = rowNames
colnames(neightborMatrixRevised) = colNames

upper_tri <- get_upper_tri(neightborMatrixRevised)
melted_Matrix <- melt(upper_tri, na.rm = TRUE)
melted_Matrix$value <- round(melted_Matrix$value,2)


plot = ggplot(data = melted_Matrix, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "black", 
                       limit = c(0,1), 
                       breaks = c(0,0.5,1),
                       labels = c(0,0.5,1),
                       space = "Lab", 
                       name="Mixing\nCoefficient") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
  theme_classic((base_size = 26)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(colour = 'black', size = 1.5),
    axis.text.x = element_text(angle = 90),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(colour = "black", size = 1.5),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.8),
    legend.direction = "horizontal",
    text=element_text(family="Helvetica"))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  scale_y_discrete(position = "right")+
  theme(legend.position = "none")
ggsave(plot, file = paste0(plot5Directory, 'FM01_mixingCoefficientMatrixTwinsAll.svg'), width =8, height = 8)

################################################################################################
###Figure 2E,F#################
################################################################################################
##################################All Big Colonies from S1######################################
################################################################################################

finalPCA1 = jointPCA %>% filter(sampleNum=="S1")
finalPCA1Big = finalPCA1 %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >=50) %>% select(-nColony) 
finalPCA1Big = finalPCA1Big %>% mutate(barcodeName = c(1:nrow(finalPCA1Big)))
finalPCA1Big$barcodeName = sub("^", "B", finalPCA1Big$barcodeName)
finalPCA1Rename = inner_join(finalPCA1,finalPCA1Big, by = "BC50StarcodeD8") %>% select(-sampleNum,-nUMI, -BC50StarcodeD8, -cellID)

neightborMatrix3 = matrix(, nrow = length(unique(finalPCA1Rename$barcodeName)), ncol = length(unique(finalPCA1Rename$barcodeName)))

for (i in c(1:length(unique(finalPCA1Rename$barcodeName)))) {
  for (j in c(1:length(unique(finalPCA1Rename$barcodeName)))) {
    barcode = unique(finalPCA1Rename$barcodeName)
    iBarcodeBx = finalPCA1Rename %>% filter(barcodeName == barcode[i])
    jBarcodeBx = finalPCA1Rename %>% filter(barcodeName == barcode[j])
    ijBarcodeBx = bind_rows(iBarcodeBx,jBarcodeBx)
    BarcodeRef = ijBarcodeBx %>% mutate(num = c(1:nrow(ijBarcodeBx)))
    S1index = BarcodeRef %>% filter(barcodeName == barcode[i]) %>% select(num)
    S2index = BarcodeRef %>% filter(barcodeName == barcode[j]) %>% select(num)
    
    knnPCA = nn2(ijBarcodeBx[,1:50], ijBarcodeBx[,1:50])
    neighborKNN = as_tibble(knnPCA[[1]]) ####gets the neighbor
    nS1 = c();
    nS2 = c();
    fraction2 = (nrow(S2index)/nrow(S1index))
    for (k in c(1:nrow(ijBarcodeBx))) {
      nS1[k] =  (sum(neighborKNN[k,2:pmin(10,nrow(S2index))] %in% S1index$num))*(nrow(S2index)/nrow(S1index))  #this part is changed
      nS2[k] =  sum(neighborKNN[k,2:pmin(10,nrow(S2index))] %in% S2index$num)
    }
    neighbor = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx)))
    neighbor = inner_join(neighbor, BarcodeRef, by = "num")
    neighborS1 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx))) %>% filter(num %in% S1index$num)
    neighborS2 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx))) %>% filter(num %in% S2index$num)
    S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(isSelf = "withSelf")  ###if want to include individual   S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(relation = "S1-S1", isSelf = "withSelf")
    S2S2 = as_tibble(neighborS2$fractionS2) %>% mutate(isSelf = "withSelf")
    S1S2 = as_tibble((neighborS1$fractionS2)) %>% mutate(isSelf ="withOther")
    S2S1 = as_tibble((neighborS2$fractionS1)) %>% mutate(isSelf = "withOther")
    toPlot = bind_rows(S1S1,S2S2,S1S2,S2S1) %>% rename(fraction = value) %>% mutate(isSelf = factor(isSelf, levels=c("withSelf","withOther")))
    withSelf = toPlot %>% filter(isSelf == "withSelf") %>% select(-isSelf)
    withOther = toPlot %>% filter(isSelf == "withOther") %>% select(-isSelf)
    neightborMatrix3[i,j]= mean(withOther$fraction)/mean(withSelf$fraction)
  }
}
rownames(neightborMatrix3) = c("B1","B2", "B3", "B4","B5", "B6", "B7", "B8","B9","B10","B11","B12","B13")
colnames(neightborMatrix3) = c("B1","B2", "B3", "B4","B5", "B6", "B7", "B8","B9","B10","B11","B12","B13")

#rownames(neightborMatrix3) = unique(finalPCA1Rename$barcodeName)
#colnames(neightborMatrix3) = unique(finalPCA1Rename$barcodeName)

meltedMatrix3 = melt(neightborMatrix3)
ggplot(data = meltedMatrix3, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() 

upper_tri <- get_upper_tri(neightborMatrix3)
melted_Matrix <- melt(upper_tri, na.rm = TRUE)
melted_Matrix$value <- round(melted_Matrix$value,2)

plot = ggplot(data = melted_Matrix, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "black", 
                       limit = c(0,1), 
                       breaks = c(0,0.5,1),
                       labels = c(0,0.5,1),
                       space = "Lab", 
                       name="Mixing\nCoefficient") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
  theme_classic((base_size = 26)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(colour = 'black', size = 1.5),
    axis.text.x = element_text(angle = 90),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(colour = "black", size = 1.5),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.8),
    legend.direction = "horizontal",
    text=element_text(family="Helvetica"))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  scale_y_discrete(position = "right")+
  theme(legend.position = "none")
ggsave(plot, file = paste0(plot5Directory, 'FM01_mixingCoefficientMatrixBigColoniesBar_S1_V1.svg'), width =8, height = 8)

#########################################################################################################
#########################################################################################################
##################################All Small/Medium Colonies from S1######################################
#########################################################################################################
finalPCA1Mid= finalPCA1 %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony <50 & nColony >30) %>% select(-nColony) 
finalPCA1Mid = finalPCA1Mid %>% mutate(barcodeName = c(1:nrow(finalPCA1Mid)))
finalPCA1Mid$barcodeName = sub("^", "B", finalPCA1Mid$barcodeName)
finalPCA1Rename = inner_join(finalPCA1,finalPCA1Mid, by = "BC50StarcodeD8") %>% select(-sampleNum,-nUMI, -BC50StarcodeD8, -cellID)

neightborMatrix4 = matrix(, nrow = length(unique(finalPCA1Rename$barcodeName)), ncol = length(unique(finalPCA1Rename$barcodeName)))
for (i in c(1:length(unique(finalPCA1Rename$barcodeName)))) {
  for (j in c(1:length(unique(finalPCA1Rename$barcodeName)))) {
    barcode = unique(finalPCA1Rename$barcodeName)
    iBarcodeBx = finalPCA1Rename %>% filter(barcodeName == barcode[i])
    jBarcodeBx = finalPCA1Rename %>% filter(barcodeName == barcode[j])
    ijBarcodeBx = bind_rows(iBarcodeBx,jBarcodeBx)
    BarcodeRef = ijBarcodeBx %>% mutate(num = c(1:nrow(ijBarcodeBx)))
    S1index = BarcodeRef %>% filter(barcodeName == barcode[i]) %>% select(num)
    S2index = BarcodeRef %>% filter(barcodeName == barcode[j]) %>% select(num)
    
    knnPCA = nn2(ijBarcodeBx[,1:50], ijBarcodeBx[,1:50])
    neighborKNN = as_tibble(knnPCA[[1]]) ####gets the neighbor
    nS1 = c();
    nS2 = c();
    fraction2 = (nrow(S2index)/nrow(S1index))
    for (k in c(1:nrow(ijBarcodeBx))) {
      nS1[k] =  (sum(neighborKNN[k,2:pmin(10,nrow(S2index))] %in% S1index$num))*(nrow(S2index)/nrow(S1index))  #this part is changed
      nS2[k] =  sum(neighborKNN[k,2:pmin(10,nrow(S2index))] %in% S2index$num)
    }
    neighbor = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx)))
    neighbor = inner_join(neighbor, BarcodeRef, by = "num")
    neighborS1 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx))) %>% filter(num %in% S1index$num)
    neighborS2 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx))) %>% filter(num %in% S2index$num)
    S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(isSelf = "withSelf")  ###if want to include individual   S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(relation = "S1-S1", isSelf = "withSelf")
    S2S2 = as_tibble(neighborS2$fractionS2) %>% mutate(isSelf = "withSelf")
    S1S2 = as_tibble((neighborS1$fractionS2)) %>% mutate(isSelf ="withOther")
    S2S1 = as_tibble((neighborS2$fractionS1)) %>% mutate(isSelf = "withOther")
    toPlot = bind_rows(S1S1,S2S2,S1S2,S2S1) %>% rename(fraction = value) %>% mutate(isSelf = factor(isSelf, levels=c("withSelf","withOther")))
    withSelf = toPlot %>% filter(isSelf == "withSelf") %>% select(-isSelf)
    withOther = toPlot %>% filter(isSelf == "withOther") %>% select(-isSelf)
    neightborMatrix4[i,j]= mean(withOther$fraction)/mean(withSelf$fraction)
  }
}
rowNames = c(1:length(unique(finalPCA1Rename$barcodeName)))
rowNames = sub("^", "B", rowNames)
colNames = rowNames
rownames(neightborMatrix4) = rowNames
colnames(neightborMatrix4) = colNames

#rownames(neightborMatrix4) = unique(finalPCA1Rename$barcodeName)
#colnames(neightborMatrix4) = unique(finalPCA1Rename$barcodeName)


upper_tri <- get_upper_tri(neightborMatrix4)
melted_Matrix <- melt(upper_tri, na.rm = TRUE)
melted_Matrix$value <- round(melted_Matrix$value,2)

plot = ggplot(data = melted_Matrix, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "black", 
                       limit = c(0,1), 
                       breaks = c(0,0.5,1),
                       labels = c(0,0.5,1),
                       space = "Lab", 
                       name="Mixing\nCoefficient") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
  theme_classic((base_size = 26)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(colour = 'black', size = 1.5),
    axis.text.x = element_text(angle = 90),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(colour = "black", size = 1.5),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.8),
    legend.direction = "horizontal",
    text=element_text(family="Helvetica"))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  scale_y_discrete(position = "right")+
  theme(legend.position = "none")
ggsave(plot, file = paste0(plot5Directory, 'FM01_mixingCoefficientMatrixMediumColoniesBar_S1_V1.svg'), width =8, height = 8)



###############################################################################################################################
###############################################################################################################################
##################################################### Sample 2 sanity check ##################################################

finalPCA1 = jointPCA %>% filter(sampleNum=="S2")
#finalPCA1Big = finalPCA1 %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony <50 & nColony >35) %>% select(-nColony) for 20-30  ####for medium colony analysis
finalPCA1Big = finalPCA1 %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >=50) %>% select(-nColony) 
finalPCA1Big = finalPCA1Big %>% mutate(barcodeName = c(1:nrow(finalPCA1Big)))
finalPCA1Big$barcodeName = sub("^", "B", finalPCA1Big$barcodeName)
finalPCA1Rename = inner_join(finalPCA1,finalPCA1Big, by = "BC50StarcodeD8") %>% select(-sampleNum,-nUMI, -BC50StarcodeD8, -cellID)

neightborMatrix5 = matrix(, nrow = length(unique(finalPCA1Rename$barcodeName)), ncol = length(unique(finalPCA1Rename$barcodeName)))

for (i in c(1:length(unique(finalPCA1Rename$barcodeName)))) {
  for (j in c(1:length(unique(finalPCA1Rename$barcodeName)))) {
    barcode = unique(finalPCA1Rename$barcodeName)
    iBarcodeBx = finalPCA1Rename %>% filter(barcodeName == barcode[i])
    jBarcodeBx = finalPCA1Rename %>% filter(barcodeName == barcode[j])
    ijBarcodeBx = bind_rows(iBarcodeBx,jBarcodeBx)
    BarcodeRef = ijBarcodeBx %>% mutate(num = c(1:nrow(ijBarcodeBx)))
    S1index = BarcodeRef %>% filter(barcodeName == barcode[i]) %>% select(num)
    S2index = BarcodeRef %>% filter(barcodeName == barcode[j]) %>% select(num)
    
    knnPCA = nn2(ijBarcodeBx[,1:50], ijBarcodeBx[,1:50])
    neighborKNN = as_tibble(knnPCA[[1]]) ####gets the neighbor
    nS1 = c();
    nS2 = c();
    fraction2 = (nrow(S2index)/nrow(S1index))
    for (k in c(1:nrow(ijBarcodeBx))) {
      nS1[k] =  (sum(neighborKNN[k,2:pmin(10,nrow(S2index))] %in% S1index$num))*(nrow(S2index)/nrow(S1index))  #this part is changed
      nS2[k] =  sum(neighborKNN[k,2:pmin(10,nrow(S2index))] %in% S2index$num)
    }
    neighbor = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx)))
    neighbor = inner_join(neighbor, BarcodeRef, by = "num")
    neighborS1 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx))) %>% filter(num %in% S1index$num)
    neighborS2 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx))) %>% filter(num %in% S2index$num)
    S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(isSelf = "withSelf")  ###if want to include individual   S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(relation = "S1-S1", isSelf = "withSelf")
    S2S2 = as_tibble(neighborS2$fractionS2) %>% mutate(isSelf = "withSelf")
    S1S2 = as_tibble((neighborS1$fractionS2)) %>% mutate(isSelf ="withOther")
    S2S1 = as_tibble((neighborS2$fractionS1)) %>% mutate(isSelf = "withOther")
    toPlot = bind_rows(S1S1,S2S2,S1S2,S2S1) %>% rename(fraction = value) %>% mutate(isSelf = factor(isSelf, levels=c("withSelf","withOther")))
    withSelf = toPlot %>% filter(isSelf == "withSelf") %>% select(-isSelf)
    withOther = toPlot %>% filter(isSelf == "withOther") %>% select(-isSelf)
    neightborMatrix5[i,j]= mean(withOther$fraction)/mean(withSelf$fraction)
  }
}

rowNames = c(1:length(unique(finalPCA1Rename$barcodeName)))
rowNames = sub("^", "B", rowNames)
colNames = rowNames
rownames(neightborMatrix5) = rowNames
colnames(neightborMatrix5) = colNames

upper_tri <- get_upper_tri(neightborMatrix5)
melted_Matrix <- melt(upper_tri, na.rm = TRUE)
melted_Matrix$value <- round(melted_Matrix$value,2)

plot = ggplot(data = melted_Matrix, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "black", 
                       limit = c(0,1), 
                       breaks = c(0,0.5,1),
                       labels = c(0,0.5,1),
                       space = "Lab", 
                       name="Mixing\nCoefficient") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
  theme_classic((base_size = 26)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(colour = 'black', size = 1.5),
    axis.text.x = element_text(angle = 90),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(colour = "black", size = 1.5),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.8),
    legend.direction = "horizontal",
    text=element_text(family="Helvetica"))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  scale_y_discrete(position = "right")+
  theme(legend.position = "none")
ggsave(plot, file = paste0(plot5Directory, 'FM01_mixingCoefficientMatrixBigColoniesBar_S2_V1.svg'), width =8, height = 8)
#########################################################################################################
#########################################################################################################
##################################All Small/Medium Colonies from S1######################################
#########################################################################################################
finalPCA1 = jointPCA %>% filter(sampleNum=="S2")
finalPCA1Mid= finalPCA1 %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony <50 & nColony >30) %>% select(-nColony) 
finalPCA1Mid = finalPCA1Mid %>% mutate(barcodeName = c(1:nrow(finalPCA1Mid)))
finalPCA1Mid$barcodeName = sub("^", "B", finalPCA1Mid$barcodeName)
finalPCA1Rename = inner_join(finalPCA1,finalPCA1Mid, by = "BC50StarcodeD8") %>% select(-sampleNum,-nUMI, -BC50StarcodeD8, -cellID)

neightborMatrix6 = matrix(, nrow = length(unique(finalPCA1Rename$barcodeName)), ncol = length(unique(finalPCA1Rename$barcodeName)))
for (i in c(1:length(unique(finalPCA1Rename$barcodeName)))) {
  for (j in c(1:length(unique(finalPCA1Rename$barcodeName)))) {
    barcode = unique(finalPCA1Rename$barcodeName)
    iBarcodeBx = finalPCA1Rename %>% filter(barcodeName == barcode[i])
    jBarcodeBx = finalPCA1Rename %>% filter(barcodeName == barcode[j])
    ijBarcodeBx = bind_rows(iBarcodeBx,jBarcodeBx)
    BarcodeRef = ijBarcodeBx %>% mutate(num = c(1:nrow(ijBarcodeBx)))
    S1index = BarcodeRef %>% filter(barcodeName == barcode[i]) %>% select(num)
    S2index = BarcodeRef %>% filter(barcodeName == barcode[j]) %>% select(num)
    
    knnPCA = nn2(ijBarcodeBx[,1:50], ijBarcodeBx[,1:50])
    neighborKNN = as_tibble(knnPCA[[1]]) ####gets the neighbor
    nS1 = c();
    nS2 = c();
    fraction2 = (nrow(S2index)/nrow(S1index))
    for (k in c(1:nrow(ijBarcodeBx))) {
      nS1[k] =  (sum(neighborKNN[k,2:pmin(10,nrow(S2index))] %in% S1index$num))*(nrow(S2index)/nrow(S1index))  #this part is changed
      nS2[k] =  sum(neighborKNN[k,2:pmin(10,nrow(S2index))] %in% S2index$num)
    }
    neighbor = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx)))
    neighbor = inner_join(neighbor, BarcodeRef, by = "num")
    neighborS1 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx))) %>% filter(num %in% S1index$num)
    neighborS2 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx))) %>% filter(num %in% S2index$num)
    S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(isSelf = "withSelf")  ###if want to include individual   S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(relation = "S1-S1", isSelf = "withSelf")
    S2S2 = as_tibble(neighborS2$fractionS2) %>% mutate(isSelf = "withSelf")
    S1S2 = as_tibble((neighborS1$fractionS2)) %>% mutate(isSelf ="withOther")
    S2S1 = as_tibble((neighborS2$fractionS1)) %>% mutate(isSelf = "withOther")
    toPlot = bind_rows(S1S1,S2S2,S1S2,S2S1) %>% rename(fraction = value) %>% mutate(isSelf = factor(isSelf, levels=c("withSelf","withOther")))
    withSelf = toPlot %>% filter(isSelf == "withSelf") %>% select(-isSelf)
    withOther = toPlot %>% filter(isSelf == "withOther") %>% select(-isSelf)
    neightborMatrix6[i,j]= mean(withOther$fraction)/mean(withSelf$fraction)
  }
}
rowNames = c(1:length(unique(finalPCA1Rename$barcodeName)))
rowNames = sub("^", "B", rowNames)
colNames = rowNames
rownames(neightborMatrix6) = rowNames
colnames(neightborMatrix6) = colNames

#rownames(neightborMatrix4) = unique(finalPCA1Rename$barcodeName)
#colnames(neightborMatrix4) = unique(finalPCA1Rename$barcodeName)


upper_tri <- get_upper_tri(neightborMatrix6)
melted_Matrix <- melt(upper_tri, na.rm = TRUE)
melted_Matrix$value <- round(melted_Matrix$value,2)

plot = ggplot(data = melted_Matrix, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "black", 
                       limit = c(0,1), 
                       breaks = c(0,0.5,1),
                       labels = c(0,0.5,1),
                       space = "Lab", 
                       name="Mixing\nCoefficient") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
  theme_classic((base_size = 26)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(colour = 'black', size = 1.5),
    axis.text.x = element_text(angle = 90),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(colour = "black", size = 1.5),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.8),
    legend.direction = "horizontal",
    text=element_text(family="Helvetica"))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  scale_y_discrete(position = "right")+
  theme(legend.position = "none")
ggsave(plot, file = paste0(plot5Directory, 'FM01_mixingCoefficientMatrixMediumColoniesBar_S2_V1.svg'), width =8, height = 8)

#####Mean value Big: 0.31
#####Mean Value Medium: 0.37
##### Mean small (20-30): 0.273

###############################################################################################################################
###############################################################################################################################
#######################################################Overall Analysis Plot ##################################################

upperTriNonSiblingsS1Large = get_upper_tri(neightborMatrix3)
diag(upperTriNonSiblingsS1Large) = NA
upperTriNonSiblingsS1Medium = get_upper_tri(neightborMatrix4)
diag(upperTriNonSiblingsS1Medium) = NA

upperTriNonSiblingsS2Large = get_upper_tri(neightborMatrix5)
diag(upperTriNonSiblingsS2Large) = NA
upperTriNonSiblingsS2Medium = get_upper_tri(neightborMatrix6)
diag(upperTriNonSiblingsS2Medium) = NA

####For siblings
upperTriNonSiblingsTwins = get_upper_tri(neightborMatrixRevised)
diag(upperTriNonSiblingsTwins) = NA


#####Calculating values for All Siblings ########

siblingsS1S2 = tibble(mixingCoefficient = numeric()) %>% add_row(mixingCoefficient = (diag(neightborMatrix))[1:12]) %>% mutate(type = "Siblings", size = "all")
nonSiblingsS1Large = tibble(mixingCoefficient = numeric()) %>% add_row(mixingCoefficient = (melt(upperTriNonSiblingsS1Large,na.rm = TRUE))$value) %>% mutate(type = "nonSiblingsS1", size = "large")
nonSiblingsS1Medium = tibble(mixingCoefficient = numeric()) %>% add_row(mixingCoefficient = (melt(upperTriNonSiblingsS1Medium,na.rm = TRUE))$value) %>% mutate(type = "nonSiblingsS1", size = "medium")
nonSiblingsS2Large = tibble(mixingCoefficient = numeric()) %>% add_row(mixingCoefficient = (melt(upperTriNonSiblingsS2Large,na.rm = TRUE))$value) %>% mutate(type = "nonSiblingsS2", size = "large")
nonSiblingsS2Medium = tibble(mixingCoefficient = numeric()) %>% add_row(mixingCoefficient = (melt(upperTriNonSiblingsS2Medium,na.rm = TRUE))$value) %>% mutate(type = "nonSiblingsS2", size = "medium")
nonSiblingsTwin = tibble(mixingCoefficient = numeric()) %>% add_row(mixingCoefficient = (melt(upperTriNonSiblingsTwins,na.rm = TRUE))$value) %>% mutate(type = "nonSiblingsTwins", size = "all")

allSiblingLargeMedium = bind_rows(siblingsS1S2, nonSiblingsS1Large, nonSiblingsS1Medium,nonSiblingsS2Large, nonSiblingsS2Medium,nonSiblingsTwin)


plot2 = ggplot(allSiblingLargeMedium, aes(x = type, y=mixingCoefficient, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  scale_fill_manual(values=c("grey", "grey", "grey", "goldenrod3")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none")
ggsave(plot2, file = paste0(plot5Directory, 'FM01_OverallPlot_V1.svg'), width = 5, height = 8)
ggsave(plot2, file = paste0(plot5Directory, 'FM01_OverallPlot_V2.svg'), width = 5, height = 8)

###############################################################################################################################
################################################UMAP PLOTS Figure 2D#####################################################################
###############################################################################################################################

umapCoordinates = as_tibble(read.table(file = paste0(home1Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum=="S1")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum=="S2")

jointPCA = inner_join(linCountTooverlaps, umapCoordinates, by = c("cellID", "sampleNum")) %>% select(-nLineages)
#BarcodesRename = jointPCA %>% select(BC50StarcodeD8) %>% unique() 
##################################################################################################
jointPCA1Barcodes = jointPCA %>% filter(sampleNum=="S1") %>% select(BC50StarcodeD8) %>% unique()
jointPCA2Barcodes = jointPCA %>% filter(sampleNum=="S2")  %>% select(BC50StarcodeD8) %>% unique()

jointBarcodesOnlyBoth = inner_join(jointPCA2Barcodes,jointPCA1Barcodes, by = "BC50StarcodeD8") 
jointBarcodesOnlyBoth = jointBarcodesOnlyBoth %>% mutate(barcodeName = c(1:nrow(jointBarcodesOnlyBoth)))
jointBarcodesOnlyBoth$barcodeName = sub("^", "B", jointBarcodesOnlyBoth$barcodeName)

jointBarcodesOnlyBothList = jointBarcodesOnlyBoth$BC50StarcodeD8
jointPCAOnlyBoth = jointPCA %>% filter(BC50StarcodeD8 %in% jointBarcodesOnlyBothList)
finalPCAJoint = inner_join(jointPCAOnlyBoth,jointBarcodesOnlyBoth, by = "BC50StarcodeD8") %>% select(-cellID, -BC50StarcodeD8, -nUMI)  ####table of PCs with SampleNum (1), BarcodeName (52)

finalPCAS1 = finalPCAJoint %>% filter(sampleNum=="S1")
finalPCAS2 = finalPCAJoint %>% filter(sampleNum=="S2")

finalPCAS1Big = finalPCAS1 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony >2) %>% select(-nColony)
finalPCAS2Big = finalPCAS2 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony >2) %>% select(-nColony)
commonS1S2Big = inner_join(finalPCAS1Big, finalPCAS2Big)
commonS1S2Big = commonS1S2Big$barcodeName
finalPCAJointBig = finalPCAJoint %>% filter(barcodeName %in% commonS1S2Big)

######using the plot from this path to decide the siblings to show
##/Users/yogesh/Dropbox (RajLab)/FateMap/data/10X/10X/2019_FM01/Analysis/plots/seurat/memory/s1s2/Delaunay/DelaunaySetColonies.pdf

finalUMAPJointtoPlot1 = finalPCAJointBig %>% filter(barcodeName == "B2")
finalUMAPJointtoPlot2 = finalPCAJointBig %>% filter(barcodeName == "B8")

plot<- ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  geom_point(data = finalUMAPJointtoPlot1, aes(x = UMAP_1, y = UMAP_2, color = sampleNum), size = 2.5, shape = 16) +
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_sibling2.svg'), width = 6, height = 6)

plot<- ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  geom_point(data = finalUMAPJointtoPlot2, aes(x = UMAP_1, y = UMAP_2, color = sampleNum), size = 2.5, shape = 16) +
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_sibling8.svg'), width = 6, height = 6)

finalPCAS1Small = finalPCAS1 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony <3) %>% select(-nColony)
finalPCAS2Small = finalPCAS2 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony <3) %>% select(-nColony)
commonS1S2small = inner_join(finalPCAS1Small, finalPCAS2Small)
commonS1S2small = commonS1S2small$barcodeName
finalPCAJointSmall = finalPCAJoint %>% filter(barcodeName %in% commonS1S2small)
#### to check:
###ggplot() + geom_point(data= jointUMAP, aes(UMAP_1, UMAP_2), color="grey1") + geom_point(data= finalPCAJointSmall, aes(UMAP_1, UMAP_2, color = sampleNum)) + facet_wrap(~barcodeName, ncol = 5) +theme_classic()

finalUMAPJointtoPlot1 = finalPCAJointSmall %>% filter(barcodeName == "B24")
finalUMAPJointtoPlot2 = finalPCAJointSmall %>% filter(barcodeName == "B64")

####size in illustrator: 96 by 84.23

plot<- ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  geom_point(data = finalUMAPJointtoPlot1, aes(x = UMAP_1, y = UMAP_2, color = sampleNum), size = 2.5, shape = 16) +
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_sibling24.svg'), width = 6, height = 6)

plot<- ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  geom_point(data = finalUMAPJointtoPlot2, aes(x = UMAP_1, y = UMAP_2, color = sampleNum), size = 2.5, shape = 16) +
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_sibling64.svg'), width = 6, height = 6)


plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_sample1UMAP.svg'), width = 6, height = 6)

plot<- ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_sample1and2UMAP.svg'), width = 6, height = 6)

plot<- ggplot() +
  geom_point(data = umapCoordinatesS2, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_sample2UMAP.svg'), width = 6, height = 6)

###############################################################################################################################
################################################UMAP PLOTS Figure 2E Non Siblings#####################################################################
###############################################################################################################################
umapCoordinates = as_tibble(read.table(file = paste0(home1Directory, "umapCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
linCountTooverlapsS1 = linCountTooverlaps %>% filter(sampleNum=="S1")
umapCoordinatesS1 = umapCoordinates %>% filter(sampleNum=="S1")
umapCoordinatesS2 = umapCoordinates %>% filter(sampleNum=="S2")

jointUMAP = inner_join(linCountTooverlaps, umapCoordinatesS1, by = c("cellID", "sampleNum")) %>% select(-nLineages)

coloniesToAnalyzeBig = jointUMAP %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >=100)
coloniesToAnalyzeMedium = jointUMAP %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony <50 & nColony >40)
coloniesToAnalyzeSmall = jointUMAP %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony <11 & nColony >2)

UMAPcoloniesToAnalyzeBig = inner_join(jointUMAP,coloniesToAnalyzeBig, by = "BC50StarcodeD8")
UMAPcoloniesToAnalyzeMedium = inner_join(jointUMAP,coloniesToAnalyzeMedium, by = "BC50StarcodeD8")
UMAPcoloniesToAnalyzeSmall = inner_join(jointUMAP,coloniesToAnalyzeSmall, by = "BC50StarcodeD8")

plot<- ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  geom_point(data = UMAPcoloniesToAnalyzeBig, aes(x = UMAP_1, y = UMAP_2, color = BC50StarcodeD8), size = 2.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "none",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_coloniesPaintedBig100.svg'), width = 6, height = 5.348)

plot<- ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  geom_point(data = UMAPcoloniesToAnalyzeMedium, aes(x = UMAP_1, y = UMAP_2, color = BC50StarcodeD8), size = 2.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "none",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_coloniesPaintedmedium4050.svg'), width = 6, height = 5.348)

plot<- ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  geom_point(data = UMAPcoloniesToAnalyzeSmall, aes(x = UMAP_1, y = UMAP_2, color = BC50StarcodeD8), size = 2.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "none",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_coloniesPaintedmedium1102.svg'), width = 6, height = 5.348)


###############################################################################################################################
################################################UMAP PLOTS Figure 2E Non Siblings#####################################################################
###############################################################################################################################
pairA = jointUMAP %>% filter( BC50StarcodeD8 %in% c("ATTCGAGTTCCTGTTGTAGGTCTTGCAGTAGTTCCAGGACTAGCACAAGG", "ATTCGTGGACGACTTCAACTTCAACGTGATGGACGAGCAGTTGCTGATGT"))
pairB = jointUMAP %>% filter( BC50StarcodeD8 %in% c("ATTCGAGTTCCTGTTGTAGGTCTTGCAGTAGTTCCAGGACTAGCACAAGG", "ATTCTAGTTGTAGTACTAGATGATCATCATGTTGTTCTTGGTCTTGTTCA"))
pairC = jointUMAP %>% filter( BC50StarcodeD8 %in% c("ATTCGAGTTCCTGTTGTAGGTCTTGCAGTAGTTCCAGGACTAGCACAAGG", "ATTCCTCCTGGTGTAGCACTACGTGCAGCTGTACATGTTGATCAAGTTGT"))  

plot<- ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  geom_point(data = pairA, aes(x = UMAP_1, y = UMAP_2, color = BC50StarcodeD8), size = 2.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "none",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_pairA.svg'), width = 6, height = 5.346)

plot<- ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  geom_point(data = pairB, aes(x = UMAP_1, y = UMAP_2, color = BC50StarcodeD8), size = 2.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "none",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_pairB.svg'), width = 6, height = 5.346)

plot<- ggplot() +
  geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93") +
  geom_point(data = pairC, aes(x = UMAP_1, y = UMAP_2, color = BC50StarcodeD8), size = 2.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "none",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(plot, file = paste0(plot5Directory, 'FM01_pairC.svg'), width = 6, height = 5.346)

#####NonSiblings

B8 = "ATTCCAGCACGACTACGTGCTGCACTTGGACCTGGTGGAGATCTACCTGG"
B6 = "ATACTTCATCAACAACCTCCTCTTCAAGTAGTACATGGTGCACTAGCAGA"
B10 = "ATTCGACGAGAACATGTAGCTGTTGTTCAACTTGGTCAAGTTCATCATGC"
B2 = "ATACGTCGTGAACTTCATCTTGGTGTTCTAGTTCCTCTACTACGACATGG"
B5 = "ATACGTGCACTTCCTGTTCCAGCTGCTGCAGAACATGGAGTAGTTCATCT"
B18 = "ATTGTACAAGATGAACAAGTACATGCACTTGGAGGTGGTCTTGTACGACA"


###########################################################################
#########################Just to plot relevant UMAPs########################
jointUMAP = inner_join(linCountTooverlaps, umapCoordinates, by = c("cellID", "sampleNum")) %>% select(-nLineages)
finalUMAPS1 = jointUMAP %>% filter(sampleNum=="S1")
finalUMAPS1Subset = finalUMAPS1 %>% filter(BC50StarcodeD8 %in% c(B6,B10))
finalUMAPS1Subset = finalUMAPS1 %>% filter(BC50StarcodeD8 %in% c(B8,B10))


