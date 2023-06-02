library(ggplot2)
library(tidyverse)
library(reshape)
library(svglite)

dataDirectory <- "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/extractedData/mouseTumor/"
plotDirectory <- "/Users/yogesh/Dropbox (RajLab)/FateMap/paper/plotData/analysisPlots/"
TumorData = as_tibble(read.table(file = paste0(dataDirectory, "tumorSummary_Exp1.csv"), header = TRUE, stringsAsFactors=F, sep = ","))

TumorDataMean = TumorData %>% group_by(group) %>% summarise(across(everything(), mean))
TumorDataSEM = TumorData %>% group_by(group) %>% summarise(across(everything(), var)) 
TumorDataSEM[,2:11] = sqrt((TumorDataSEM[,2:11]/(nrow(TumorData)/2))*(((nrow(TumorData)/2)-1))/(nrow(TumorData)/2)) ##population variance

TumorDataMean = melt(as.data.frame(TumorDataMean)) %>% dplyr::rename(mean = value)
TumorDataMean$variable = sub("X", "", TumorDataMean$variable)
TumorDataSEM = melt(as.data.frame(TumorDataSEM))  %>% dplyr::rename(sem = value)
TumorDataSEM$variable = sub("X", "", TumorDataSEM$variable)

TumorDataAll = as_tibble(inner_join(TumorDataMean,TumorDataSEM, by = c("group", "variable")))
TumorDataAll$variable = as.numeric(TumorDataAll$variable)

plot <- ggplot(TumorDataAll, aes(x=variable, y=mean, group=group)) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.y = element_blank(), 
        axis.text.x=element_blank())

ggsave(plot = plot, file = paste0(plotDirectory, 'TumorDataAllDose','.svg'), width = 6, height = 4)

#####t-test 
tTest = tibble()

pValue2 = tibble(pvalue = numeric())
  
for (i in 2:ncol(TumorData)) {
  tTest = t.test(TumorData[11:20,i], TumorData[1:10,i], alternative = "greater",conf.level = 0.9)
  pValue2[i-1,1] = tTest[["p.value"]]
}

write.table(pValue2, file=paste0(plotDirectory,'pValueMouseTumorFig4.tsv'), col.names = TRUE, sep='\t')

#######
# Latest EDIT:YG 2/20/2022 for revisions
data2Directory <- "/Users/ygy1258/Dropbox (RajLab)/FateMap/paper/extractedData/mouseTumor/"
plot2Directory <- "/Users/ygy1258/Dropbox (RajLab)/FateMap/paper/plotData/revisionAnalysisPlots/"

TumorData = as_tibble(read.table(file = paste0(data2Directory, "tumorSummary_Exp2.csv"), header = TRUE, stringsAsFactors=F, sep = ","))

TumorDataMean = TumorData %>% group_by(group) %>% summarise(across(everything(), mean))
TumorDataSEM = TumorData %>% group_by(group) %>% summarise(across(everything(), var))
TumorDataSEM[,2:12] = sqrt((TumorDataSEM[,2:12]/(nrow(TumorData)/2))*(((nrow(TumorData)/2)-1))/(nrow(TumorData)/2)) ##population variance

TumorDataMean = melt(as.data.frame(TumorDataMean)) %>% dplyr::rename(mean = value)
TumorDataMean$variable = sub("X", "", TumorDataMean$variable)
TumorDataSEM = melt(as.data.frame(TumorDataSEM))  %>% dplyr::rename(sem = value)
TumorDataSEM$variable = sub("X", "", TumorDataSEM$variable)

TumorDataAll = as_tibble(inner_join(TumorDataMean,TumorDataSEM, by = c("group", "variable")))
TumorDataAll$variable = as.numeric(TumorDataAll$variable)

plot <- ggplot(TumorDataAll, aes(x=variable, y=mean, group=group)) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  theme_classic((base_size = 24)) + theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = "black", size = 2), text=element_text(family="Helvetica")) +
  theme(axis.title.y = element_blank(), 
        axis.text.x=element_blank())
ggsave(plot = plot, file = paste0(plot2Directory, 'TumorDataAllDose','.svg'), width = 8, height = 4)

tTest = tibble()

pValue2 = tibble(pvalue = numeric())

for (i in 2:ncol(TumorData)) {
  tTest = t.test(TumorData[1:5,i], TumorData[6:10,i], alternative = "greater",conf.level = 0.9)
  pValue2[i-1,1] = tTest[["p.value"]]
}

write.table(pValue2, file=paste0(plot2Directory,'revision_pValueMouseTumor.tsv'), col.names = TRUE, sep='\t')

