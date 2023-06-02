library(tidyverse)
library(msm)
library(reshape2)
library(tictoc)

###directories
plotDirectory <- "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/plotData/analysisPlots/"


####calculating the frequency of efficiency

 
  
#initializing parameters
cellNumber = 200000
barcodeNumber = 50000000
MOI = 0.12 #average rate of success NOT GFP+% (based on experimental data)
divisionNumberCol = c(4,5)
splitNumber = 4
cellcultureLoss = 0.05
injectionLoss = 0
tumorextractionLoss = 0
dnaextractionLoss = 0
libraryprepLoss = .05
repetitions = 1000

totalBarcodedCells =  20000+20000 # from transduction charts
survivingCellsTotal = 87+14+14 ###from BC02 result
reprogEffic = survivingCellsTotal/totalBarcodedCells

barcodes <- c(1:barcodeNumber)


results <- list()

for (i in 1:repetitions) {
  #SEEDING
  tic(paste("repetition", i))
  
  initialCells <-
    data.frame(
      "Cell ID" = 1:cellNumber,
      "Barcode" = rep(0, cellNumber),
      "Daughters" = rep(0, cellNumber),
      "Split ID" = rep(0, cellNumber),
      "Loss1" = rep(FALSE, cellNumber),
      "Loss2" = rep(FALSE, cellNumber),
      "Reprogram" = rep(FALSE, cellNumber)
    )
  
  #BARCODING
  sampledBarcodes <- sample(barcodes, cellNumber, replace = TRUE)
  
  poissonCutoff <- 1 - dpois(0, MOI, log = FALSE)
  barcodeIndex <- runif(cellNumber, 0, 1) < poissonCutoff
  
  sampledBarcodes[barcodeIndex == FALSE] <- 0
  initialCells$Barcode <- sampledBarcodes
  
  #VARIABLE DIVISIONS
  iterativeCells <- initialCells
  
  divisionNumber = sample(divisionNumberCol,1)
  for (j in 1:divisionNumber) {
    sampledDaughters <-
      round(rtnorm(nrow(iterativeCells), 2, 0.25, 0))  ###what is rtnorm? Why truncated normal distribution? and how were mean and sd decided?
    # sampledDaughters <- round(rtnorm(nrow(iterativeCells), 1.75, 0.5, 0))
    iterativeCells$Daughters <- sampledDaughters
    
    iterativeCells <-
      iterativeCells[rep(row.names(iterativeCells), iterativeCells$Daughters), 1:7]
  }
  
  #SPLITTING
  splits <- c(1:splitNumber)
  sampledSplits <-
    sample(splits, nrow(iterativeCells), replace = TRUE)
  
  iterativeCells$Split.ID <- sampledSplits
  
  #PREDRUG LOSS
  loss_predrug = cellcultureLoss + injectionLoss
  lossIndex1 <- runif(nrow(iterativeCells), 0, 1) < loss_predrug
  iterativeCells$Loss1 <- lossIndex1
  
  finalCells <- filter(iterativeCells, Loss1 == FALSE)
  
  #PLATE SPLITTING
  plate1 <- filter(finalCells, Split.ID == 1)
  plate2 <- filter(finalCells, Split.ID == 2)
  
  #REPROGRAMMING
  
  plate1Max <- round(nrow(plate1) * reprogEffic)
  plate2Max <- round(nrow(plate2) * reprogEffic)
  
  reprogIndex1 <- runif(nrow(plate1), 0, 1) < reprogEffic
  plate1$Reprogram <- reprogIndex1
  reprogIndex2 <- runif(nrow(plate2), 0, 1) < reprogEffic
  plate2$Reprogram <- reprogIndex2
  
  plate1ReprogSampled <- filter(plate1, Reprogram == TRUE)
  plate2ReprogSampled <- filter(plate2, Reprogram == TRUE)
  
  #POSTDRUG LOSS
  loss_postdrug = tumorextractionLoss + dnaextractionLoss + libraryprepLoss
  
  lossIndex21 <-
    runif(nrow(plate1ReprogSampled), 0, 1) < loss_postdrug
  plate1ReprogSampled$Loss2 <- lossIndex21
  
  lossIndex22 <-
    runif(nrow(plate2ReprogSampled), 0, 1) < loss_postdrug
  plate2ReprogSampled$Loss2 <- lossIndex22
  
  plate1 <- filter(plate1, Loss2 == FALSE)
  plate2 <- filter(plate2, Loss2 == FALSE)
  
  #OVERLAP ANALYSIS
  plate1ReprogBarcodes <-
    unique(select(filter(plate1, Barcode != 0, Reprogram == TRUE),
                  Barcode))
  plate2ReprogBarcodes <-
    unique(select(filter(plate2, Barcode != 0, Reprogram == TRUE),
                  Barcode))
  
  if (nrow(union(plate1ReprogBarcodes, plate2ReprogBarcodes)) == 0) {
    overlap <- 0
  } else{
    overlap <-
      nrow(intersect(plate1ReprogBarcodes, plate2ReprogBarcodes)) / nrow(union(plate1ReprogBarcodes, plate2ReprogBarcodes))
  }
  
  results <- append(results, overlap)
  
  toc()
}

resultsMelt <- melt(results)
colnames(resultsMelt) <- c("Percent.Overlap","Heritability")
write.table(resultsMelt, file=paste0(plotDirectory,'randomSimulationInVitro.tsv'), col.names = TRUE, sep='\t')

expt1Overlap = 100*87/(87+14+14)



plot = ggplot(resultsMelt*100, aes(x=Percent.Overlap)) + geom_histogram(binwidth=1) +theme_classic() + 
  geom_vline(aes(xintercept=expt1Overlap),color="blue", size=1) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 2), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica"))

ggsave(plot, file = paste0(plotDirectory, 'BC02_MemoryRandomExperimentalDataInVitro.svg'), width = 4, height = 2)
