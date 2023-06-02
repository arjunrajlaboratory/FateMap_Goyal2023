These tables were made using the biomaRt package in R, using the code below. The "hg19" table used the ENSG gene_ids from the hg19 gtf file used by the pipeline as of June 2019. The code below was run on June 18th, 2019 to create an internal database of Ensembl gene_id to HGNC gene_symbol IDs. The hg38 table was created on July 30th, 2019. To make a new one for a new reference genome, you can edit this block of R code: (note that this code is also in a comment block in the normalizeMeltedCountsAndAddGeneSymbol.R script)


library(biomaRt)
library(tidyverse)
library(GenomicFeatures)
library(here)
gtfFileUsedByPipeline <- '/Users/emsanford/Downloads/Homo_sapiens.GRCh38.97.gtf'
txdb <- makeTxDbFromGFF(file = gtfFileUsedByPipeline, format="gtf")
vectorOfENSG_IDs_to_find_HGNC_symbols_for <- names(sum(width(IRanges::reduce(exonsBy(txdb2, by = "gene")))))  ## this is just a vector of "ENSG0000XXXX" strings. This line of code gives you this vector from a parsed GTF file that you make below
outputTableWithGeneNameSymbols <- here('newGeneNameSymbolTable.tsv') 
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "uswest.ensembl.org", ensemblRedirect = FALSE)
hgnc.name.table <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
                         filters = 'ensembl_gene_id', 
                         values = vectorOfENSG_IDs_to_find_HGNC_symbols_for, 
                         mart = ensembl)

hgnc.name.tibble <- tibble(ensg = hgnc.name.table$ensembl_gene_id, hgnc_symbol = hgnc.name.table$hgnc_symbol)
write_tsv(hgnc.name.tibble, outputTableWithGeneNameSymbols, col_names = TRUE)