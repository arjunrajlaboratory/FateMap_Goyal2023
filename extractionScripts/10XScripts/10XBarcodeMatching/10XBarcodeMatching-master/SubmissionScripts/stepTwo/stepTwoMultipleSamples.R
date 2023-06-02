#======================================================================================================================================
#Input Files: It takes inputs from files generated from stepTwo.R for each sample --> stepTwoCellIDUMIBarcodes.txt
#Specify the PATH containing folder of individual samples in HomeDirectory. 
#Sub folders containing samples should be specified in sampleFolders
#======================================================================================================================================

homeDirectory <- '/project/arjunrajlab/experimentName/Analysis/stepTwo/'
sampleFolders = c('s1/', 's2/', 's3/', 's4/')

#***************************************************************************************************************************************
#*****************************************************DO NOT EDIT BEYOND THIS POINT*****************************************************
#***************************************************************************************************************************************
library(dplyr)
library(tidyverse)

numberSamples = length(sampleFolders);
stepTwoCellIDUMIBarcodesAll = list()


for(i in 1:numberSamples){
  stepTwoCellIDUMIBarcodes = as_tibble(read.table(paste0(homeDirectory, sampleFolders[i], 'stepTwoCellIDUMIBarcodes.txt'), stringsAsFactors=F, header = TRUE)) %>% mutate(sampleNum = i)
  
  if(is.null(dim(stepTwoCellIDUMIBarcodesAll))){
    stepTwoCellIDUMIBarcodesAll = stepTwoCellIDUMIBarcodes
  } else {
    stepTwoCellIDUMIBarcodesAll = bind_rows(stepTwoCellIDUMIBarcodesAll, stepTwoCellIDUMIBarcodes)
  }
  
}

write.table(stepTwoCellIDUMIBarcodesAll, file= paste0(homeDirectory,'stepTwoCellIDUMIBarcodes.txt'),row.names=F,col.names=T,quote=F,sep="\t")
write.table(stepTwoCellIDUMIBarcodesAll[,4], file= paste0(homeDirectory,'stepTwoBarcodes50.txt'),row.names=F,col.names=F,quote=F,sep="\t")
write.table(stepTwoCellIDUMIBarcodesAll[,5], file= paste0(homeDirectory,'stepTwoBarcodes40.txt'),row.names=F,col.names=F,quote=F,sep="\t")
write.table(stepTwoCellIDUMIBarcodesAll[,6], file= paste0(homeDirectory,'stepTwoBarcodes30.txt'),row.names=F,col.names=F,quote=F,sep="\t")
