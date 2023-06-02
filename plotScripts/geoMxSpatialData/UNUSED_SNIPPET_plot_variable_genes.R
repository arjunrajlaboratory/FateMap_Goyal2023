cpm_summary =cpm_long %>% filter(region == 'FullROI') %>% group_by(patient_id, gene, surg_immune_pre) %>% summarize(mean_cpm = mean(cpm), sd_cpm = sd(cpm), cv_cpm = sd(cpm) / mean(cpm))
foo=cpm_summary %>% group_by(patient_id, surg_immune_pre) %>% slice_max(cv_cpm, n = 5)

library(ggrepel)
ggplot(cpm_summary) +
  aes(x = mean_cpm, y = sd_cpm, label = gene) +
  geom_point() +
  facet_grid(patient_id ~ surg_immune_pre, scales = 'free')

ggplot(cpm_summary) +
  aes(x = mean_cpm, y = sd_cpm, label = gene) +
  geom_point() +
  facet_grid(patient_id ~ ., scales = 'free') +
  geom_text_repel(data = foo)

ggplot(cpm_long %>%
         dplyr::filter(patient_id == 163,
                       region == 'FullROI',
         )) +
  aes(x= cpm, y= cpm, color = as.factor(S100B_level), label = ROI) +
  geom_point(shape = 16) +
  facet_grid(factor(gene, levels = genes_of_interest$gene)~patient_plug,
             scales = 'free')+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
  )

ggplot(cpm_long %>%
         dplyr::filter(patient_id == 129,
                       gene %in% genes_of_interest$gene,
                       region == 'FullROI',
                       resistant_type %in% c('M', "T")
         )) +
  aes(x= full_sample_ordered, y= cpm, color = as.factor(S100B_level), label = ROI) +
  geom_point(shape = 16) +
  facet_grid(factor(gene, levels = genes_of_interest$gene)~patient_plug,
             scales = 'free')+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
  )
ggplot(cpm_long %>%
         dplyr::filter(
                       gene %in% c("S100B"),
                       region == 'FullROI',
                       surg_immune_pre == 'SURG'
         )) +
  aes(x= S100B_level, y= cpm, label = ROI) +
  geom_point(shape = 16) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
  ) +
  geom_smooth(method = lm)
