library(tidyverse)
library(here)
library(readxl)
library(edgeR)
library(cowplot)
library(lsr)
library(Seurat)
library(Matrix)
library(MASS)
library(nnls)

# data reading and wrangling ----------------------------------------------

#read input
log_cpm <- read_csv(file = here('Paper', 'extractedData', 'geomx_log_cpm.csv'))
FM01 = readRDS(here('Paper', 'extractedData', 'from_yogo', 's1s2ScTransform_50pcs_filter.rds'))
dataDirectory1 <- "~/Dropbox (RajLab)/REVISIONS/FateMapAnalysis_Maalavika/FM01_DEGandVarGenes/Data/"
#Read pre-existing DEGs from Maalavika
DEG_FM01 <- read.delim(paste0(dataDirectory1, "DEG_FC2_pval0.01.tsv"), sep = "\t")$gene

#Make signature matrix
#ACTA2, which marks smooth muscles, is found largely in Seurat cluster 8; IFIT2, which marks type-1 interferon signaling, is found largely in cluster 12; VCAM1, which marks vascular adhesion, is found largely in cluster 15; NGFR, which marks neural crest cells, is found largely in cluster 7; MLANA, which marks melanocytes, is found largely in clusters 0 and 3.

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
#Do TMM normalization
de_sig_norm <- calcNormFactors(de_sig, method = "TMM")
signature_mat <- log_cpm(de_sig_norm, log=T)

#START FUNCTION
# ROIs of interest- 13/14
dat_log_cpm = log_cpm[,colnames(log_cpm) %in% c('BRAF-MEKI-2_013_FullROI', 'BRAF-MEKI-2_014_FullROI')]
dat_log_cpm = dat_log_cpm %>% as_tibble(.) %>% rowwise() %>% mutate(roi_mean = mean(c(`BRAF-MEKI-2_013_FullROI`, `BRAF-MEKI-2_014_FullROI`), na.rm = T))
dat_log_cpm$gene = rownames(log_cpm)

# Assuming bulk_data is a matrix with samples in columns and genes in rows
# Assuming signature_matrix is a matrix with cell types in columns and genes in rows
# Both matrices should have row names containing the gene names

dat_log_cpm_mat = as.matrix(dat_log_cpm[,1:3])
rownames(dat_log_cpm_mat) = dat_log_cpm[['gene']]

# Align the genes between the bulk_data and the signature_matrix
common_genes <- intersect(rownames(dat_log_cpm_mat), rownames(signature_mat))
common_genes = intersect(common_genes, DEG_FM01)
dat_log_cpm_mat <- dat_log_cpm_mat[common_genes,]
dat_log_cpm_mat = na.omit(dat_log_cpm_mat)
signature_mat <- signature_mat[common_genes,]
signature_mat = na.omit(signature_mat)


# Perform NNLS deconvolution
deconvoluted_data <- matrix(0, nrow = ncol(dat_log_cpm_mat), ncol = ncol(signature_mat))
for (i in 1:ncol(dat_log_cpm_mat)) {
  nnls_result <- nnls(A = signature_mat, b = dat_log_cpm_mat[, i])
  deconvoluted_data[i,] <- nnls_result$x
}

# Normalize the proportions to sum up to 1
cell_proportions <- t(apply(deconvoluted_data, 1, function(x) x / sum(x)))

number_nuclei = metadata %>% filter(full_sample %in% colnames(dat_log_cpm_mat)) %>% pull(AOINucleiCount)
number_nuclei[3] = sum(number_nuclei)
cell_numbers = matrix(0, nrow = nrow(cell_proportions), ncol = ncol(cell_proportions))
for (i in 1:nrow(cell_proportions)) {
  cell_numbers[i,] = cell_proportions[i,] * number_nuclei[i]
}

cell_numbers = round(cell_numbers)
pseudobulk_data2 = tibble(pseudobulk_roi_1 = rowSums(sweep(signature_mat, 2, as.numeric(cell_numbers[1,]), "*"))/sum(cell_numbers[1,]),
                          pseudobulk_roi_2 = rowSums(sweep(signature_mat, 2, as.numeric(cell_numbers[2,]), "*"))/sum(cell_numbers[2,]))

#original, true distance
true_distance = dist(t(dat_log_cpm_mat[,1:2]))

#reconstructed distance, to use as 'observed' distance
observed_distance = dist(t(pseudobulk_data2))

do_sampling = function(x){
  count_vector = cell_numbers[3,]
  cell_pool = rep(1:length(count_vector), times = count_vector)
  sampled_idx = sample(x = cell_pool, size = number_nuclei[1], replace = F)
  # Count the occurrences of each sampled index
  sampled_indices_tibble <- as_tibble(table(sampled_idx), .name_repair = 'minimal')
  sampled_indices_tibble = sampled_indices_tibble %>% rename(cell_type = sampled_idx, roi_1 = n)
  sampled_indices_tibble = as_tibble(count_vector) %>% mutate(cell_type = as.character(1:length(count_vector))) %>% full_join(sampled_indices_tibble, by = join_by(cell_type)) %>% rename(total_cells = value)
  sampled_indices_tibble = sampled_indices_tibble %>% mutate(roi_2 = total_cells - roi_1)
  sampled_indices_tibble[is.na(sampled_indices_tibble)] = 0
  
  pseudobulk_data = tibble(pseudobulk_roi_1 = rowSums(sweep(signature_mat, 2, sampled_indices_tibble$roi_1, "*"))/sum(sampled_indices_tibble$roi_1),
         pseudobulk_roi_2 = rowSums(sweep(signature_mat, 2, sampled_indices_tibble$roi_2, "*"))/sum(sampled_indices_tibble$roi_2))
  
  #Distance metric, consider changing
  cur_distance = dist(rbind(pseudobulk_data[[1]], pseudobulk_data[[2]]))
  #cur_distance = cosine(pseudobulk_data[[1]], pseudobulk_data[[2]])
  return(cur_distance)
}
set.seed(42)
distances = sapply(1:10000, do_sampling)
hist(distances)

combined_cell_proportion = cell_numbers[3,]

roi1_cells = expand.grid(seq(from = 0, to = combined_cell_proportion[1], by = 10),
         seq(from = 0, to = combined_cell_proportion[2], by = 10),
         seq(from = 0, to = combined_cell_proportion[3], by = 10),
         seq(from = 0, to = combined_cell_proportion[4], by = 10),
         seq(from = 0, to = combined_cell_proportion[5], by = 10))
roi2_cells = -1*sweep(roi1_cells, 2, combined_cell_proportion, FUN = "-")

#Restrict to samples that are within 2% of the original cell numbers
roi1_good_idx = rowSums(roi1_cells) <= number_nuclei[1] + 0.02 * number_nuclei[1] & rowSums(roi1_cells) >= number_nuclei[1] - 0.02 * number_nuclei[1]
roi1_cells = roi1_cells[roi1_good_idx,]
roi2_cells = roi2_cells[roi1_good_idx,]

unequal_proportion_distances = 0

for (i in 1:nrow(roi1_cells)){
  pseudo1 = rowSums(sweep(signature_mat, 2, as.numeric(roi1_cells[i,]),"*"))/sum(roi1_cells[i,])
  pseudo2 = rowSums(sweep(signature_mat, 2, as.numeric(roi2_cells[i,]),"*"))/sum(roi2_cells[i,])
  unequal_proportion_distances[i] = dist(rbind(pseudo1, pseudo2))
}

plot_data = tibble(roi1_cells) 
plot_data = sweep(plot_data, 2, combined_cells, '/')
plot_data = plot_data %>% select_if(~ !any(is.na(.)))


plot_data = plot_data %>%
  mutate(total_abs_deviation = rowSums(abs(. - 0.5)),
         mad = rowMeans(abs(. - 0.5)))

plot_data = plot_data %>% add_column(distance = unequal_proportion_distances)

#END FUNCTION
ggplot(plot_data) +
  aes(x= total_abs_deviation, y = distance) +
  geom_point(alpha= 0.3, size = 0.3) +
  geom_smooth()+
  xlab('cell proportion euclidean distance') +
  ylab('tx distance')

ggplot(plot_data) +
  aes(x= total_abs_deviation, fill = distance > observed_distance) +
  geom_histogram(position = 'identity', alpha = 0.5) + 
  xlab('cell proportion euclidean distance') +
  ylab('tx distance')

ggplot(plot_data) +
  aes(x= mad, fill = distance > observed_distance) +
  geom_histogram(position = 'identity', alpha = 0.5) + 
  xlab('cell proportion euclidean distance') +
  ylab('tx distance')
