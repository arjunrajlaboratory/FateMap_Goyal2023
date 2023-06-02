
#Load libraries that will be used in the scripts.
library(tidyverse)
library(gplots)
library(RColorBrewer)


rm(list=ls())

##########
#Combining two melted datasets together#
##########

#set the working directory. 
setwd("/Users/ygy1258/Dropbox (RajLab)/FateMap/paper/extractedData/bulkRNASeq/")

#deRandomizing#
##########

RNASeqCounts = as_tibble(read.table(file = "BulkColonySeq_20200220_meltedDataAll_normalized_hg19.tsv", header = TRUE, sep = "\t"))
Random_demultiplex_Revised = as_tibble(read.table(file = "RandomizationRevised.tsv", header = TRUE, sep = "\t"))

##
total = RNASeqCounts %>% dplyr::select(sampleID,counts) %>% group_by(sampleID) %>% summarise(counts = sum(counts)) %>% arrange(counts)
total
##
Random_demultiplex_Revised$Original = sub("^", "sample_", Random_demultiplex_Revised$Original)

RNASeqCounts1 = RNASeqCounts %>% mutate(sampleID = gsub("WM989-ResistantColonies-","",RNASeqCounts$sampleID),
                                        sampleID = as.integer(sampleID),
                                        rawCorrespondence = sub("^", "WM989_RC", sampleID)) 

RNASeqCounts1$rawCorrespondence[RNASeqCounts1$rawCorrespondence == "WM989_RC1"] <-  "WM989_RC01"
RNASeqCounts1$rawCorrespondence[RNASeqCounts1$rawCorrespondence == "WM989_RC2"] <-  "WM989_RC02"
RNASeqCounts1$rawCorrespondence[RNASeqCounts1$rawCorrespondence == "WM989_RC3"] <-  "WM989_RC03"
RNASeqCounts1$rawCorrespondence[RNASeqCounts1$rawCorrespondence == "WM989_RC4"] <-  "WM989_RC04"
RNASeqCounts1$rawCorrespondence[RNASeqCounts1$rawCorrespondence == "WM989_RC5"] <-  "WM989_RC05"
RNASeqCounts1$rawCorrespondence[RNASeqCounts1$rawCorrespondence == "WM989_RC6"] <-  "WM989_RC06"
RNASeqCounts1$rawCorrespondence[RNASeqCounts1$rawCorrespondence == "WM989_RC7"] <-  "WM989_RC07"
RNASeqCounts1$rawCorrespondence[RNASeqCounts1$rawCorrespondence == "WM989_RC8"] <-  "WM989_RC08"
RNASeqCounts1$rawCorrespondence[RNASeqCounts1$rawCorrespondence == "WM989_RC9"] <-  "WM989_RC09"


RNASeqCountsRevised = inner_join(RNASeqCounts1,Random_demultiplex_Revised, by = "sampleID") %>% dplyr::select(-sampleID) %>% arrange(Original)
write_tsv(RNASeqCountsRevised, "/Users/ygy1258/Dropbox (RajLab)/FateMap/paper/extractedData/bulkRNASeq/deRandomized_BulkColonySeq_20200220_meltedData_hg19_normalizedRevised.tsv", col_names = T)

