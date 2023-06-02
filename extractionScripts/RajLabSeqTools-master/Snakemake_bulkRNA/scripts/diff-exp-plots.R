### Differential-expression-plotting

# Last updated 2019-08-29, P. Burnham

# Will produce two plots to compare bulk RNA seq data.
# Volcano plot based on the two conditions.
# Heatmaps with significance based on two groups

### Input
args = commandArgs(trailingOnly=TRUE) # arg.1 = tsv file name, arg.2 = # genes for heatmap,
                                      # arg.3 = log10.significance.threshold, arg.4 = log2-fold-threshold.

# extract info from argument 1

base.name = gsub(pattern = ".diffexp.tsv", replacement = "", x = rev(unlist(strsplit(args[1],split = "/")))[1]) ;
base.path = paste(rev(rev(unlist(strsplit(args[1],split = "/")))[-1]),collapse = "/") ;
condition.base = unlist(strsplit(base.name,split = "-"))[3] ;
condition.altered = unlist(strsplit(base.name,split = "-"))[1] ;

# Generate Ffgures
source("scripts/differential-expression-custom-functions.R")

results.table = read.table(args[1], header = T) ;

# heatmap for expression
pdf(file=paste0(base.path,"/",base.name,".diffexp.heatmap.pdf"),  width=6, height=6,  pointsize=12)

Heatmap_DEseq(results.table = results.table, selected.gene.number = as.numeric(args[2]), 
              condition.base = condition.base, condition.altered = condition.altered) 

dev.off()

# volcano plot
pdf(file=paste0(base.path,"/",base.name,".diffexp.volcano.pdf"),  width=6,  height=4, pointsize=12)

Volcano_DEseq(results.table = results.table, significance = 10**as.numeric(args[3]), 
              log2.fold.cutoff = as.numeric(args[4]), label.pts = T) 

dev.off()