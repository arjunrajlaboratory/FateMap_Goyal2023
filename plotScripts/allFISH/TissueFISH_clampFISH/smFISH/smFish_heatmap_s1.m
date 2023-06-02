% Heatmap of smFISH Data (FateMap) %%
%This script creates heatmap from the smFISH data in YG102s1.

clear

% import gene count table for s1
formatSpec = '%f%f%f%f%f';
spotTable= readtable('FateMap/extractedData/smFISH/YG102s1/geneCounts2_s1.csv','Format',formatSpec);
% check
head(spotTable,5)


%% make heatmap for each gene

% ACTA 2
figure(1)
heatmap_gene1 = heatmap(spotTable,"X","Y","ColorVariable","ACTA2", "ColorMethod","none")
heatmap_gene1.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene1,'FateMap/extractedData/smFISH/YG102s1/ACTA2_heatmap2_s1.tif')
% BGN
figure(2)
heatmap_gene2 = heatmap(spotTable,"X","Y","ColorVariable","BGN", "ColorMethod","none")
heatmap_gene2.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene2,'FateMap/extractedData/smFISH/YG102s1/BGN_heatmap2_s1.tif')
% UNC
figure(3)
heatmap_gene3 = heatmap(spotTable,"X","Y","ColorVariable","UBC", "ColorMethod","none")
heatmap_gene3.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene3,'FateMap/extractedData/smFISH/YG102s1/UBC_heatmap2_s1.tif')

%% normalize spot counts to housekeeping gene UBC **optional**

spotTable.normACTA2 = spotTable.ACTA2./spotTable.UBC;
spotTable.normBGN = spotTable.BGN./spotTable.UBC;
spotTable.normUBC = spotTable.UBC./spotTable.UBC;
% check
spotTable(1:500:end,:)

%% make normalized heatmap for each gene

% ACTA 2
figure(4)
heatmap_gene1 = heatmap(spotTable,"X","Y","ColorVariable","normACTA2", "ColorMethod","none")
heatmap_gene1.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene1,'FateMap/extractedData/smFISH/YG102s1/normACTA2_heatmap2_s1.tif')
% BGN
figure(5)
heatmap_gene2 = heatmap(spotTable,"X","Y","ColorVariable","normBGN", "ColorMethod","none")
heatmap_gene2.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene2,'FateMap/extractedData/smFISH/YG102s1/normBGN_heatmap2_s1.tif')
% UNC
figure(6)
heatmap_gene3 = heatmap(spotTable,"X","Y","ColorVariable","normUBC", "ColorMethod","none")
heatmap_gene3.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene3,'FateMap/extractedData/smFISH/YG102s1/normUBC_heatmap2_s1.tif')
