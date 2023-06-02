library(tidyverse)
library(here)
library(edgeR)
library(lsr)
library(Seurat)
library(Matrix)
library(MASS)
library(nnls)
library(cowplot)
source(here('Paper', 'deconvolution_function.R'))
# Read input
log_cpm <- read.csv(file = here('Paper', 'extractedData', 'geomx_log_cpm.csv'), row.names = 1)
FM01 = readRDS(here('Paper', 'extractedData', 'from_yogo', 's1s2ScTransform_50pcs_filter.rds'))
metadata = read_csv(file = here('Paper', 'extractedData', 'processed_metadata.csv'))
# data reading and wrangling ----------------------------------------------

# Read pre-existing DEGs from Maalavika
dataDirectory1 <- "~/Dropbox (RajLab)/REVISIONS/FateMapAnalysis_Maalavika/FM01_DEGandVarGenes/Data/"
DEG_FM01 <- read.delim(paste0(dataDirectory1, "DEG_FC2_pval0.01.tsv"), sep = "\t")$gene

# Make signature matrix
FM01_sub = subset(FM01, integrated_snn_res.0.6 %in% c(8, 12, 15, 7, 0, 3))

case_resistant_type = function(cluster){
  case_when(
    cluster == 8 ~ 'ACTA2_type',
    cluster == 12 ~ 'IFIT2_type',
    cluster == 15 ~ 'VCAM1_type',
    cluster == 7 ~ 'NGFR_type',
    cluster %in% c(0, 3) ~ 'MLANA_type'
  )
}
FM01_sub@meta.data = FM01_sub@meta.data %>% mutate(resistant_type = case_resistant_type(integrated_snn_res.0.6))
Idents(FM01_sub) = 'resistant_type'

avg_exprs = AverageExpression(FM01_sub, return.seurat = F, assays = 'SCT', slot = 'counts')

signature_rna = avg_exprs$SCT

de_sig = DGEList(signature_rna)
# Do TMM normalization
de_sig_norm <- calcNormFactors(de_sig, method = "TMM")
signature_matrix <- cpm(de_sig_norm, log=T)

#Usage of the function
pt_134 = deconvolve_data(log_cpm, signature_matrix, metadata, c('BRAF-MEKI-2_013_FullROI', 'BRAF-MEKI-2_014_FullROI'))
pt_163 = deconvolve_data(log_cpm, signature_matrix, metadata, c('BRAF-MEKI-2_076_FullROI', 'BRAF-MEKI-2_077_FullROI'))
pt_129 = deconvolve_data(log_cpm, signature_matrix, metadata, c('BRAF-MEKI-2_045_FullROI', 'BRAF-MEKI-2_046_FullROI'))
pt_109 = deconvolve_data(log_cpm, signature_matrix, metadata, c('BRAF-MEKI-2_025_FullROI', 'BRAF-MEKI-2_026_FullROI'))


h1=ggplot() +
  aes(x = pt_134$distances) +
  geom_histogram() +
  geom_vline(xintercept = pt_134$observed_distance) +
  labs(subtitle = paste0('P(shuffled distance >= observed distance) = ', pt_134$null_pvalue))
h2=ggplot() +
  aes(x = pt_163$distances) +
  geom_histogram() +
  geom_vline(xintercept = pt_163$observed_distance) +
  labs(subtitle = paste0('P(shuffled distance >= observed distance) = ', pt_163$null_pvalue))
h3=ggplot() +
  aes(x = pt_129$distances) +
  geom_histogram() +
  geom_vline(xintercept = pt_129$observed_distance) +
  labs(subtitle = paste0('P(shuffled distance >= observed distance) = ', pt_129$null_pvalue))
h4=ggplot() +
  aes(x = pt_109$distances) +
  geom_histogram() +
  geom_vline(xintercept = pt_109$observed_distance) +
  labs(subtitle = paste0('P(shuffled distance >= observed distance) = ', pt_109$null_pvalue))

th1 = plot_grid(h1,h2,h3,h4, nrow =4)
ggsave(here('Paper', 'plots', 'null_histograms.svg'), th1)
StatPercentileX <- ggproto("StatPercentileX", Stat,
                           compute_group = function(data, scales, probs) {
                             percentiles <- quantile(data$x, probs=probs)
                             data.frame(xintercept=percentiles)
                           },
                           required_aes = c("x")
)

stat_percentile_x <- function(mapping = NULL, data = NULL, geom = "vline",
                              position = "identity", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    stat = StatPercentileX, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

StatPercentileXLabels <- ggproto("StatPercentileXLabels", Stat,
                                 compute_group = function(data, scales, probs) {
                                   percentiles <- quantile(data$x, probs=probs)
                                   data.frame(x=percentiles, y=Inf,
                                              label=paste0("p", probs*100, ": ",
                                                           round(percentiles, digits=3)))
                                 },
                                 required_aes = c("x")
)

stat_percentile_xlab <- function(mapping = NULL, data = NULL, geom = "text",
                                 position = "identity", na.rm = FALSE,
                                 show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    stat = StatPercentileXLabels, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

p1=ggplot(pt_134$plot_data) +
  aes(x= mad, fill = distance > pt_134$observed_distance) +
  geom_histogram(position = 'identity', alpha = 0.5) + 
  xlab('mean fraction away from equal partition for each cell type') +
  ylab('tx distance') +
  stat_percentile_x(probs=c(0.25), linetype=2) +
  stat_percentile_xlab(probs=c(0.25), hjust=1, vjust=1.5, angle=90)

p2=ggplot(pt_109$plot_data) +
  aes(x= mad, fill = distance > pt_109$observed_distance) +
  geom_histogram(position = 'identity', alpha = 0.5) + 
  xlab('mean fraction away from equal partition for each cell type') +
  ylab('tx distance')+
  stat_percentile_x(probs=c(0.25), linetype=2) +
  stat_percentile_xlab(probs=c(0.25), hjust=1, vjust=1.5, angle=90)

p3=ggplot(pt_129$plot_data) +
  aes(x= mad, fill = distance > pt_129$observed_distance) +
  geom_histogram(position = 'identity', alpha = 0.5) + 
  xlab('mean fraction away from equal partition for each cell type') +
  ylab('tx distance')+
  stat_percentile_x(probs=c(0.25), linetype=2) +
  stat_percentile_xlab(probs=c(0.25), hjust=1, vjust=1.5, angle=90)

p4=ggplot(pt_163$plot_data) +
  aes(x= mad, fill = distance > pt_163$observed_distance) +
  geom_histogram(position = 'identity', alpha = 0.5) + 
  xlab('mean fraction away from equal partition for each cell type') +
  ylab('tx distance')+
  stat_percentile_x(probs=c(0.25), linetype=2) +
  stat_percentile_xlab(probs=c(0.25), hjust=1, vjust=1.5, angle=90)

tp1 = plot_grid(p1, p2, p3, p4, nrow=4, align = 'hv')
ggsave(here('Paper', 'plots', 'unequal_proportion_histograms.svg'), tp1)
