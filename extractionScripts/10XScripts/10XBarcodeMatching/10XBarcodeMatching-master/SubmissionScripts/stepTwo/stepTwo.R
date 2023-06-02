#======================================================================================================================================
#Input Files: It takes inputs from files generated in stepOne and from 10XCellranger pipeline generated filteredMatrix -> barcode.tsv.gz
#Change these file PATHs based on your folder structure and where your datasets are stored. 
#$PATH needs to be changed/input at 3 places.
#======================================================================================================================================

input1Directory <- '/project/arjunrajlab/10Xdatasets/experimentName/s1/outs/filtered_feature_bc_matrix/'
input2Directory <- '/project/arjunrajlab/experimentName/Analysis/stepOne/s1/'
outputDirectory <- '/project/arjunrajlab/experimentName/Analysis/stepTwo/s1/'

#***************************************************************************************************************************************
#*****************************************************DO NOT EDIT BEYOND THIS POINT*****************************************************
#***************************************************************************************************************************************

library(stringdist)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyverse)

data1file = as_tibble(read.table(paste0(input1Directory,"barcodes.tsv.gz"), stringsAsFactors=F)) %>% dplyr::rename(cellID = V1) 
data1file = as_tibble(substring(data1file$cellID, 1,nchar(data1file[1,1])-2)) %>% dplyr::rename(cellID = value) 
data2file = as_tibble(read.table(paste0(input2Directory,"uniqueShavedReads.txt"), stringsAsFactors=F)) %>% dplyr::rename(cellID = V1,UMI = V2, BC = V3) %>%
  mutate(BC50 = substring(BC,1,50),
         BC40 = substring(BC,1,40),
         BC30a = substring(BC,1,30),
         BC30b = substring(BC,1,30))
cellIDUMIBarcodes = inner_join(data1file, data2file, by = "cellID")
Barcodes = unique(cellIDUMIBarcodes$BC)
cellIDs = unique(cellIDUMIBarcodes$cellID)

set.seed(2059)
subsample1 = sample(Barcodes,5000)
subsample2 = sample(Barcodes,5000)
subsample3 = sample(Barcodes,5000)
BarcodesLv1 = as.integer(stringdistmatrix(subsample1, method = "lv"))
BarcodesLv2 = as.integer(stringdistmatrix(subsample2, method = "lv"))
BarcodesLv3 = as.integer(stringdistmatrix(subsample3, method = "lv"))
lBarcodesLv = length(BarcodesLv1)

BarcodesLv = tibble(
  lvdist = c(BarcodesLv1, BarcodesLv2, BarcodesLv3),
  subsamNum = c(rep("subsamping1", lBarcodesLv), rep("subsamping2", lBarcodesLv), rep("subsamping3", lBarcodesLv)))

BarcodesLvHist <- BarcodesLv %>% group_by(subsamNum, lvdist) %>% summarise(length(lvdist)) %>%
  group_by(subsamNum) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

BarcodesLvHistPlot <- ggplot(BarcodesLvHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  facet_wrap(facets = vars(subsamNum)) +
  theme_classic()

set.seed(2059)
cellIDsLv = tibble(lvdist = as.integer(stringdistmatrix(cellIDs, method = "lv")))
cellIDsHist <- cellIDsLv  %>% group_by(lvdist)%>% summarise(length(lvdist)) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

cellIDsHistPlot <- ggplot(cellIDsHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  theme_classic()

#writing files
ggsave(BarcodesLvHistPlot,file=paste0(outputDirectory,'stepTwoBarcodesLvBeforeStarcode.pdf'))
ggsave(cellIDsHistPlot,file=paste0(outputDirectory,'stepTwoCellIdsLvBeforeStarcode.pdf'))
write.table(cellIDUMIBarcodes, file= paste0(outputDirectory,'stepTwoCellIDUMIBarcodes.txt'),row.names=F,col.names=T,quote=F,sep="\t")
write.table(cellIDUMIBarcodes[,4], file= paste0(outputDirectory,'stepTwoBarcodes50.txt'),row.names=F,col.names=F,quote=F,sep="\t")
write.table(cellIDUMIBarcodes[,5], file= paste0(outputDirectory,'stepTwoBarcodes40.txt'),row.names=F,col.names=F,quote=F,sep="\t")
write.table(cellIDUMIBarcodes[,6], file= paste0(outputDirectory,'stepTwoBarcodes30.txt'),row.names=F,col.names=F,quote=F,sep="\t")



