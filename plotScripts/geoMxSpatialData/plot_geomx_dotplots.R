library(tidyverse)
library(here)
library(readxl)
library(edgeR)
library(cowplot)

# data reading and wrangling ----------------------------------------------

#read input
roi_input = read_csv(here('Data', 'DSP', 'processed', 'GEOMX-ROI.csv'))
filtered_data = read_xlsx(here('Data', 'DSP', 'processed', '5_filter_RB.xlsx'), sheet = 3)
segment_prop = read_xlsx(here('Data', 'DSP', 'processed', '5_filter_RB.xlsx'), sheet = 1)
FM_deg_list = read.delim(here('Paper', 'extractedData', 'from_yogo', 'CloneGenes_FC1_p0.05_AllSinglets.tsv'), row.names = 1)
#wrangle the technical properties of each ROI
segment_prop = segment_prop %>% dplyr::select(-c(SlideName, ScanLabel, SegmentDisplayName))
segment_prop$ROILabel = as.numeric(segment_prop$ROILabel)
segment_prop$SegmentLabel = str_remove(segment_prop$SegmentLabel, ' ')
#wrangle input data
sample_names = colnames(filtered_data)[-1]
sample_names = str_remove_all(sample_names, coll(' ')) %>% str_replace_all(fixed('|'), fixed('_'))
colnames(filtered_data) = c('gene', sample_names)
#wrangle metadata
metadata = tibble('full_sample' = sample_names)
metadata = metadata %>% separate(full_sample, into = c(NA, 'ROI', 'region'), sep = '_', remove = F)
metadata$ROI = as.numeric(metadata$ROI)
metadata = metadata %>% left_join(roi_input)
metadata = metadata %>% left_join(segment_prop, by = c('ROI' = 'ROILabel', 'region' = 'SegmentLabel'))
metadata = metadata %>% unite("patient_plug", patient_letter, plug, remove = F)
write_csv(metadata, file = here('Paper', 'extractedData', 'processed_metadata.csv'))

#Make DESeq object
gene_list = filtered_data[[1]]
filtered_data_mat = as.matrix(filtered_data[-1])
rownames(filtered_data_mat) = gene_list
de = DGEList(filtered_data_mat)
#Do TMM normalization
de_norm <- calcNormFactors(de, method = "TMM")

#gene list from figure 1 of FateMap paper (plus AXL and BGN and immune genes)
genes_of_interest = tibble::tribble(
  ~resistant_type,    ~gene,
  "A",  "ACTA2",
  "A",  "ACTG2",
  "A",  "MYOCD",
  "A",  "TAGLN",
  "B",  "IFIT2",
  "B",   "OASL",
  "B",   "CCL3",
  "B",  "DDX58",
  "B",    "AXL",
  "C",  "VCAM1",
  "C",  "PKDCC",
  "C",   "TDO2",
  "C",  "FOXF2",
  "D",   "NGFR",
  "D", "COL9A3",
  "D",  "S100B",
  "D",  "ITGA6",
  "D",    "BGN",
  "E",  "MLANA",
  "E",  "SOX10",
  "E",   "MITF",
  "E",   "PMEL",
  "M",   "CD14",
  "M",   "CD68",
  "T",    "CD2",
  "T",   "CD3E",
  "T",    "CD4",
  "T",    "CD8A"
)

limited_genes_of_interest = tibble::tribble(
  ~resistant_type,    ~gene,
  "1",    "AXL",
  "A",  "ACTA2",
  "B",  "IFIT2",
  "B",   "OASL",
  "C",  "VCAM1",
  "D",   "NGFR",
  "E",  "MLANA",
  "E",  "SOX10"
)

# make long cpm for dotplots
cpm <- cpm(de_norm, log=F)
write.csv(cpm, file = here('Paper', 'extractedData', 'geomx_cpm.csv'), row.names = T, col.names = T)
cpm_log <- cpm(de_norm, log=T)
write.csv(cpm_log, file = here('Paper', 'extractedData', 'geomx_log_cpm.csv'), row.names = T, col.names = T)

cpm_long = as_tibble(t(cpm), rownames = 'full_sample')
cpm_long = cpm_long %>% pivot_longer(cols = -full_sample, names_to = 'gene', values_to = 'cpm')
cpm_long = cpm_long %>% left_join(metadata)
cpm_long = cpm_long %>% left_join(genes_of_interest)
cpm_long$full_sample_ordered = reorder(cpm_long$full_sample, cpm_long$S100B_level)
cpm_summary =cpm_long %>% filter(region == 'FullROI') %>% group_by(patient_id, gene, surg_immune_pre) %>% summarize(mean_cpm = mean(cpm), sd_cpm = sd(cpm), cv_cpm = sd(cpm) / mean(cpm))

write_csv(cpm_summary, file = here('Paper', 'extractedData', 'geomx_variable_genes_by_patient_condition.csv'))

# MAIN FIGURE PLOTS, two patients, 7 genes ---------------------------------------------------------------------

p1=ggplot(cpm_long %>%
         dplyr::filter(patient_id == 163,
                       patient_plug != 'NA_5',
                       gene %in% limited_genes_of_interest$gene,
                       region == 'FullROI'
         )) +
  aes(x= full_sample_ordered, y= cpm, color = as.factor(S100B_level), label = ROI) +
  geom_point(shape = 16) +
  facet_grid(factor(gene, levels = limited_genes_of_interest$gene)~patient_plug,
             scales = 'free')+
  theme_cowplot()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
  )
ggsave(here('Paper', 'plots', 'yogo_revision', 'patient_163_rev2.eps'), plot = p1, width = 10, height = 12)
# plug 4 and 5 are pre-treatment
p2=ggplot(cpm_long %>%
         dplyr::filter(patient_id == 134,
                       gene %in% limited_genes_of_interest$gene,
                       region == 'FullROI'
         )) +
  aes(x= full_sample_ordered, y= cpm, color = as.factor(S100B_level), label = ROI) +
  geom_point(shape = 16) +
  facet_grid(factor(gene, levels = limited_genes_of_interest$gene)~factor(patient_plug, levels = c('NA_4',
                                                                                           'NA_5',
                                                                                           'NA_1',
                                                                                           'NA_2',
                                                                                           'NA_3',
                                                                                           'NA_6',
                                                                                           'NA_7',
                                                                                           'NA_8',
                                                                                           'NA_9',
                                                                                           'NA_10',
                                                                                           'NA_11')),
             scales = 'free')+
  theme_cowplot()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
  )
plot_grid(p1+guides(color='none'),p2, align = 'hv', nrow = 1, rel_widths = c(1.2,3))
ggsave(here('Paper', 'plots', 'yogo_revision', 'combined.eps'), width = 9.5, height = 5.5)

# SUPPLEMENTAL FIGURE,  all patients and all genes ------------------------

s1=ggplot(cpm_long %>%
         dplyr::filter(patient_id == 163,
                       gene %in% genes_of_interest$gene,
                       region == 'FullROI'
         )) +
  aes(x= full_sample_ordered, y= cpm, color = as.factor(S100B_level), label = ROI) +
  geom_point(shape = 16) +
  facet_grid(factor(gene, levels = genes_of_interest$gene)~patient_plug,
             scales = 'free')+
  theme_cowplot()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
  )

s2=ggplot(cpm_long %>%
         dplyr::filter(patient_id == 134,
                       gene %in% genes_of_interest$gene,
                       region == 'FullROI'
         )) +
  aes(x= full_sample_ordered, y= cpm, color = as.factor(S100B_level), label = ROI) +
  geom_point(shape = 16) +
  facet_grid(factor(gene, levels = genes_of_interest$gene)~factor(patient_plug, levels = c('NA_4',
                                                                                                   'NA_5',
                                                                                                   'NA_1',
                                                                                                   'NA_2',
                                                                                                   'NA_3',
                                                                                                   'NA_6',
                                                                                                   'NA_7',
                                                                                                   'NA_8',
                                                                                                   'NA_9',
                                                                                                   'NA_10',
                                                                                                   'NA_11')),
             scales = 'free')+
  theme_cowplot()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
  )

### Patient 129 has metadata with unknown meaning
s3=ggplot(cpm_long %>%
         dplyr::filter(patient_id == 129,
                       gene %in% genes_of_interest$gene,
                       region == 'FullROI'
         )) +
  aes(x= full_sample_ordered, y= cpm, color = as.factor(S100B_level)) +
  geom_point(shape = 16) +
  facet_grid(factor(gene, levels = genes_of_interest$gene)~patient_plug,
             scales = 'free')+  
  theme_cowplot()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
  )

#Plug 3 is pre-treatment
s4=ggplot(cpm_long %>%
         dplyr::filter(patient_id == 109,
                       gene %in% genes_of_interest$gene,
                       region == 'FullROI'
         )) +
  aes(x= full_sample_ordered, y= cpm, color = as.factor(S100B_level), label = ROI) +
  geom_point(shape = 16) +
  facet_grid(factor(gene, levels = genes_of_interest$gene)~factor(patient_plug, levels = c('NA_3',
                                                                                           'NA_1',
                                                                                           'NA_2')),
             scales = 'free')+
  theme_cowplot()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
  )
plot_grid(s1+guides(color='none'),
          s2+guides(color='none'),
          s4+guides(color='none'),
          s3,
          align = 'h',
          nrow = 1,
          rel_widths = c(1.5,2,1,2.8))
ggsave(here('Paper', 'plots', 'yogo_revision', 'supplemental.eps'), width = 20, height = 18)
