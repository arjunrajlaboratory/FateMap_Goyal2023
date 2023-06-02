library(tidyverse)
library(ggpubr)
library(here)
library(Seurat)
geomx_var_genes = read_csv(here('Paper', 'extractedData', 'geomx_variable_genes_by_patient_condition.csv'))
log_cpm = read_csv(here('Paper', 'extractedData', 'geomx_variable_genes_by_patient_condition.csv'))
FM01 = readRDS(here('Paper', 'extractedData', 'from_yogo', 's1s2ScTransform_50pcs_filter.rds'))
expressed_in_989 = rownames(FM01@assays$RNA@counts)[rowSums(FM01@assays$RNA@counts) > 0]

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

all_markers = FindAllMarkers(FM01_sub, only.pos = T)
all_markers = all_markers %>% dplyr::filter(p_val_adj < 0.05)
#Pick the top 100 markers by log2fc for each resistant type

top_markers = all_markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 200)
geomx_var_genes = geomx_var_genes %>% mutate(in_FM_markers = gene %in% top_markers[['gene']])
geomx_var_genes = geomx_var_genes %>% mutate(CV_percentile = percent_rank(cv_cpm))
geomx_var_genes = geomx_var_genes %>% dplyr::filter(gene %in% expressed_in_989)

top_500_geomx = geomx_var_genes %>% group_by(patient_id) %>% slice_max(order_by = cv_cpm, n = 1000)

overlaps = top_500_geomx %>% group_by(patient_id) %>% summarize(sum(gene %in% top_markers[['gene']])/1000)

p <- ggboxplot(geomx_var_genes, x = "in_FM_markers", y = "cv_cpm",
               fill = "in_FM_markers", palette = "jco",
               facet.by = c("patient_id", 'surg_immune_pre'), short.panel.labs = T)
p + stat_compare_means(vjust = .5)

ggsave(here('Paper', 'plots', 'marker_CV_comparison.svg'), plot=p)
