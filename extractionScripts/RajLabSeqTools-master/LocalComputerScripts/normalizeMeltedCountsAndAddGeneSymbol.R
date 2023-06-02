library(tidyverse)
library(GenomicFeatures)
library(here)
library(refGenome)
# note: if you get errors from any of the above "library" lines, e.g. "there is no package called ‘here’, 
#       enter the following command in the RStudio console: install.packages(<name_of_library_that_failed_to_load>)

######### here beginneth user-defined parameters #########

meltedDataInputFile   <- '/Users/emsanford/Dropbox (RajLab)/Shared_Eric/Signal_Integration/Analysis_SI2-SI4/extractedData/SI3_5sample_rerun_meltedData.tsv'
gtfFileUsedByPipeline <- '/Users/emsanford/Dropbox (RajLab)/Shared_Eric/Signal_Integration/Analysis_SI2-SI4/refs/hg38.gtf'
gtfFileDirectory      <- '/Users/emsanford/Dropbox (RajLab)/Shared_Eric/Signal_Integration/Analysis_SI2-SI4/refs'
ENSGtoGeneSymbolTable <- '/Users/emsanford/Dropbox (Personal)/Eric/Penn/raj_lab/code/rajlabseqtools/default/LocalComputerScripts/geneSymbolConversionTables/hg38_EnsgHgncSymbolMapping.tsv'
outputFile            <- '/Users/emsanford/Dropbox (RajLab)/Shared_Eric/Signal_Integration/Analysis_SI2-SI4/extractedData/SI3_5sample_rerun_RNA-seq-pipeline-output.tsv'

######### here endeth user-defined parameters ############

# parse GTF file into a first R object, used to retrieve gene names
ens <- ensemblGenome()
setwd(gtfFileDirectory)
read.gtf(ens, "hg38.gtf")
my_gene <- getGenePositions(ens)

# parse GTF file into a second R object. Use this one to make a union model from which to calculate gene lengths
txdb <- makeTxDbFromGFF(file = gtfFileUsedByPipeline, format="gtf")
lengthsPergeneid <- sum(width(IRanges::reduce(exonsBy(txdb, by = "gene"))))
lengthtbl <- tibble(gene_id = names(lengthsPergeneid), length = lengthsPergeneid)

# read in the "melted" HTSeq count table and add the following normalized values for each sample: RPM, RPKM, and TPM
htseq.table <- read_tsv(meltedDataInputFile, col_names = T)
htseq.table.withLength <- inner_join(htseq.table, lengthtbl, by = 'gene_id') # this step also removes non_gene features from the table, e.g. "__no_feature"
htseq.table.withRPM <- htseq.table.withLength %>% 
  group_by(experiment, sampleID) %>%
  mutate(totalMappedReads = sum(counts), rpm = 1000000*counts / totalMappedReads)
htseq.table.withRPKM <- htseq.table.withRPM %>%
  mutate(rpkm = 1000*rpm/length)
# note this TPM calc method doesn't take average fragment size into account (some papers use it, some papers don't)
htseq.table.withTPM <- htseq.table.withRPKM %>%
  mutate(ctsOverLength = counts/length, denominator = sum(ctsOverLength), tpm = ctsOverLength/denominator * 1000000) %>%
  dplyr::select(-c(ctsOverLength, denominator))

# now add the gene symbol to the final table and write the output file
HGNC.symbol.table <- read_tsv(ENSGtoGeneSymbolTable, col_names = T)
HGNC.symbol.table <- HGNC.symbol.table %>% mutate(gene_id = ensg) %>% dplyr::select(-ensg)
htseq.table.with.HGNCsymbol <- plyr::join(data.frame(htseq.table.withTPM), data.frame(HGNC.symbol.table), by = 'gene_id', type="left", match="first")  #uses plyr join on a data frame instead of dplyr left_join on a tibble due to first match option
htseq.table.almostfinal <- as_tibble(htseq.table.with.HGNCsymbol)

# now add the HGNC gene symbol to the final table and write the output file
HGNC.symbol.table <- read_tsv(ENSGtoGeneSymbolTable, col_names = T)
HGNC.symbol.table <- HGNC.symbol.table %>% mutate(gene_id = ensg) %>% dplyr::select(-ensg)
htseq.table.with.HGNCsymbol <- plyr::join(data.frame(htseq.table.withTPM), data.frame(HGNC.symbol.table), by = 'gene_id', type="left", match="first")  #uses plyr join on a data frame instead of dplyr left_join on a tibble due to first match option
htseq.table.almostfinal <- as_tibble(htseq.table.with.HGNCsymbol)

# add hg38.gtf names to the genes that are missing hgnc symbols
indices.genes.missing.hgnc.symbols <- which(is.na(htseq.table.almostfinal$hgnc_symbol))
gene.ids.missing.symbols <- htseq.table.almostfinal$gene_id[indices.genes.missing.hgnc.symbols]
temp.table1 <- data.frame(tibble(gene_id = gene.ids.missing.symbols))
temp.table2 <- data.frame(tibble(gene_id = my_gene$gene_id, name = my_gene$gene_name))
temp.table <- plyr::join(temp.table1, temp.table2, by = 'gene_id', type="left", match="first") 
htseq.table.almostfinal$hgnc_symbol[indices.genes.missing.hgnc.symbols] <- temp.table$name
stopifnot(sum(is.na(htseq.table.almostfinal$hgnc_symbol)) == 0)
# rename HGNC column to "gene_name" column. Gene names are the HGNC symbol if one exists for the ensembl ID and the gene_name in the gtf file otherwise.
htseq.table.almostfinal[["gene_name"]] <- htseq.table.almostfinal$hgnc_symbol
htseq.table.almostfinal <- dplyr::select(htseq.table.almostfinal, -hgnc_symbol)

htseq.table.final <- htseq.table.almostfinal %>% dplyr::select(-length, -totalMappedReads)
write_tsv(htseq.table.final, outputFile, col_names = T)


####### if you need to make a new ENSGtoGeneSymbolTable, uncomment, edit, and run this block of code ###############
# library(biomaRt)
# library(tidyverse)
# library(GenomicFeatures)
# library(here)
# gtfFileUsedByPipeline <- '/Users/emsanford/Downloads/Homo_sapiens.GRCh38.97.gtf'
# txdb <- makeTxDbFromGFF(file = gtfFileUsedByPipeline, format="gtf")
# vectorOfENSG_IDs_to_find_HGNC_symbols_for <- names(sum(width(IRanges::reduce(exonsBy(txdb2, by = "gene")))))  ## this is just a vector of "ENSG0000XXXX" strings. This line of code gives you this vector from a parsed GTF file that you make below
# outputTableWithGeneNameSymbols <- here('newGeneNameSymbolTable.tsv') 
# ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "uswest.ensembl.org", ensemblRedirect = FALSE)
# hgnc.name.table <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
#                          filters = 'ensembl_gene_id', 
#                          values = vectorOfENSG_IDs_to_find_HGNC_symbols_for, 
#                          mart = ensembl)
# 
# hgnc.name.tibble <- tibble(ensg = hgnc.name.table$ensembl_gene_id, hgnc_symbol = hgnc.name.table$hgnc_symbol)
# write_tsv(hgnc.name.tibble, outputTableWithGeneNameSymbols, col_names = TRUE)
####################################################################################################################