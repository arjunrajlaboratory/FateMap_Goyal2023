% Heatmap of smFISH Data (FateMap) %%
%This script creates heatmap from the smFISH data in YG102s5.

clear

% import gene count table for s5
formatSpec = '%f%f%f%f%f';
spotTable= readtable('FateMap/extractedData/smFISH/YG102s5/geneCounts_s5.csv','Format',formatSpec);
% check
head(spotTable,5)

%% make heatmap for each gene

% SOX10
figure(1)
heatmap_gene1 = heatmap(spotTable,"X","Y","ColorVariable","SOX10", "ColorMethod","none")
heatmap_gene1.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene1,'FateMap/extractedData/smFISH/YG102s5/SOX10_heatmap_s5.tif')
% AXL
figure(2)
heatmap_gene2 = heatmap(spotTable,"X","Y","ColorVariable","AXL", "ColorMethod","none")
heatmap_gene2.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene2,'FateMap/extractedData/smFISH/YG102s5/AXL_heatmap_s5.tif')
% NGFR
figure(3)
heatmap_gene3 = heatmap(spotTable,"X","Y","ColorVariable","NGFR", "ColorMethod","none")
heatmap_gene3.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene3,'FateMap/extractedData/smFISH/YG102s5/NGFR_heatmap_s5.tif')

