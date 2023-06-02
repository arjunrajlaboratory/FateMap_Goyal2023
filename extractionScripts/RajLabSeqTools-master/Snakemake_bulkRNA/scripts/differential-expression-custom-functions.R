# differential-expression-custom-functions.R

# Last updated: 02019-08-29, P. Burnham

# Packages
library(biomaRt) ; library(ggplot2); library(ggrepel); library(dplyr); library(reshape2); library(ggpubr)


# Functions 
Volcano_DEseq = function(results.table, significance = 10**-5 , log2.fold.cutoff = 3, label.pts = F){
  result.dataframe = data.frame(results.table)
  result.dataframe$geneID = rownames(result.dataframe)
  
  plotted = ggplot()+
    geom_point(data = result.dataframe[(abs(result.dataframe$log2FoldChange) < log2.fold.cutoff)|
                                         (result.dataframe$pvalue > significance),],
               aes(log2FoldChange, -1*log10(pvalue)),
               col="grey")+
    geom_point(data = result.dataframe[(abs(result.dataframe$log2FoldChange) >= log2.fold.cutoff)&
                                         (result.dataframe$pvalue <= significance),],
               aes(log2FoldChange, -1*log10(pvalue)), col = "red")+
    geom_vline(xintercept = c(-1*log2.fold.cutoff,log2.fold.cutoff), linetype =2)+
    geom_hline(yintercept = -1*log10(significance), linetype =2)+
    labs(x = "log2 Fold Change", y = "-log10(p-value)")+
    theme_bw()
  
  if(label.pts){
    selected.genes = result.dataframe[(abs(result.dataframe$log2FoldChange) >= log2.fold.cutoff) & (result.dataframe$pvalue <= significance),]
    
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    hgnc_convert <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),filters = 'ensembl_gene_id', 
                          values = selected.genes$geneID, mart = ensembl)
    colnames(hgnc_convert)[1] = "geneID"
    
    selected.genes = merge(selected.genes,hgnc_convert, all.x=T)
    selected.genes[is.na(selected.genes$hgnc_symbol),]$hgnc_symbol = selected.genes[is.na(selected.genes$hgnc_symbol),]$geneID
    plotted = plotted+ geom_text_repel(data = selected.genes,size=3,   segment.size = 0.25,
                                       box.padding = unit(0.35, "lines"),
                                       point.padding = unit(0.3, "lines"),
                                       aes(log2FoldChange, -1*log10(pvalue),label = hgnc_symbol))
  }
  
  return(plotted)
}

Heatmap_DEseq = function(results.table,   selected.gene.number = 15, condition.base = "State", condition.altered = "Change"){
  result.dataframe = data.frame(results.table)
  result.dataframe$geneID = rownames(result.dataframe)
  
  result.dataframe$compareMean = (result.dataframe$baseMean)*(2**(result.dataframe$log2FoldChange))
  
  most.diff.df = result.dataframe[(rank(result.dataframe$log2FoldChange,ties.method = "random") > (nrow(result.dataframe) - selected.gene.number)) | 
                                    (rank(result.dataframe$log2FoldChange,ties.method = "random") < selected.gene.number),]
  most.diff.df$baseMean.log = log(most.diff.df$baseMean)
  most.diff.df$compareMean.log = log(most.diff.df$compareMean)
  
  select.diff.df = subset(most.diff.df, select = c("geneID", "baseMean.log", "compareMean.log"))
  colnames(select.diff.df)[2:3] = c(condition.base, condition.altered)
  melted.df = melt(select.diff.df)
  colnames(melted.df)[2:3] = c("Group", "logCount")
  
  most.diff.df = merge(most.diff.df,melted.df)
  
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  hgnc_convert <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),filters = 'ensembl_gene_id', 
                        values = most.diff.df$geneID, mart = ensembl)
  colnames(hgnc_convert)[1] = "geneID"
  
  most.diff.df = merge(most.diff.df,hgnc_convert, all.x=T)
  most.diff.df[is.na(most.diff.df$hgnc_symbol) | (most.diff.df$hgnc_symbol == "") ,]$hgnc_symbol = most.diff.df[is.na(most.diff.df$hgnc_symbol)| (most.diff.df$hgnc_symbol == ""),]$geneID
  
  most.diff.df$gene = factor(most.diff.df$hgnc_symbol, levels = unique(most.diff.df[order(most.diff.df$log2FoldChange),]$hgnc_symbol))
  
  heatmap = ggplot(most.diff.df, aes(Group,gene))+
    geom_tile(aes(fill = logCount), col ="black")+
    scale_fill_gradient(low = "purple",high = "yellow")+ labs(title = "Expression")+
    theme_bw() + 
    theme(legend.position = "left",
          axis.title.x = element_text("Helvetica", size = 10),
          axis.text = element_text("Helvetica", size = 6), 
          plot.title = element_text("Helvetica", size = 14),
          axis.title.y = element_blank())
  
  significance.plot = ggplot(data = unique(subset.data.frame(most.diff.df, select = c("padj","gene"))), aes(log10(padj),gene))+
    geom_point()+scale_x_reverse()+
    theme_bw()+labs(title = "Significance", x = "Log10(p-value, adjusted) ")+
    theme(legend.position = "", 
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text("Helvetica", size = 10),
          axis.text.x = element_text("Helvetica", size = 6), 
          plot.title = element_text("Helvetica", size = 14))
  
  
  return( ggarrange(heatmap, significance.plot, align = "v", widths = c(7,2)))
}
