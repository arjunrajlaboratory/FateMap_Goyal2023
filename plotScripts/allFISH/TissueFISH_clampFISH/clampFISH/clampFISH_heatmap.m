%% Heatmap of clampFISH Data (FateMap) %%
%This script creates heatmap from the clampFISH data in E158rep2.

clear

% import gene count table
formatSpec = '%d%d%f%f%f%f%f%f%f%f%f%f%f%f';
spotTable= readtable('subregion4_geneCount2.csv','Format',formatSpec);
% check
head(spotTable,5)



%% make heatmap for each gene

% R1_YFP_500ms_UBC
figure(1)
heatmap_gene1 = heatmap(spotTable,"X","Y","ColorVariable","R1_YFP_500ms_UBC", "ColorMethod","none")
heatmap_gene1.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene1,'R1_YFP_500ms_UBC_heatmap_E158rep2.tif')

% R1_CY3_250ms_NGFR
figure(2)
heatmap_gene2 = heatmap(spotTable,"X","Y","ColorVariable","R1_CY3_250ms_NGFR", "ColorMethod","none")
heatmap_gene2.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene2,'R1_CY3_250ms_NGFR_heatmap_E158rep2.tif')

% R1_A594_500ms_MMP1
figure(3)
heatmap_gene3 = heatmap(spotTable,"X","Y","ColorVariable","R1_A594_500ms_MMP1", "ColorMethod","none")
heatmap_gene3.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene3,'R1_A594_500ms_MMP1_heatmap_E158rep2.tif')

% R1_CY5_250ms_AXL
figure(4)
heatmap_gene4 = heatmap(spotTable,"X","Y","ColorVariable","R1_CY5_250ms_AXL", "ColorMethod","none")
heatmap_gene4.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene4,'R1_CY5_250ms_AXL_heatmap_E158rep2.tif')

% R2_YFP_500ms_UBC
figure(5)
heatmap_gene5 = heatmap(spotTable,"X","Y","ColorVariable","R2_YFP_500ms_UBC", "ColorMethod","none")
heatmap_gene5.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene5,'R2_YFP_500ms_UBC_heatmap_E158rep2.tif')

% R2_CY3_250ms_ITGA3
figure(6)
heatmap_gene6 = heatmap(spotTable,"X","Y","ColorVariable","R2_CY3_250ms_ITGA3", "ColorMethod","none")
heatmap_gene6.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene6,'R2_CY3_250ms_ITGA3_heatmap_E158rep2.tif')

% R2_CY5_500ms_EGFR
figure(7)
heatmap_gene7 = heatmap(spotTable,"X","Y","ColorVariable","R2_CY5_500ms_EGFR", "ColorMethod","none")
heatmap_gene7.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene7,'R2_CY5_500ms_EGFR_heatmap_E158rep2.tif')

% R2_A594_1000ms_FN1
figure(8)
heatmap_gene8 = heatmap(spotTable,"X","Y","ColorVariable","R2_A594_1000ms_FN1", "ColorMethod","none")
heatmap_gene8.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene8,'R2_A594_1000ms_FN1_heatmap_E158rep2.tif')

% R3_YFP_500ms_UBC
figure(9)
heatmap_gene9 = heatmap(spotTable,"X","Y","ColorVariable","R3_YFP_500ms_UBC", "ColorMethod","none")
heatmap_gene9.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene9,'R3_YFP_500ms_UBC_heatmap_E158rep2.tif')

% R3_CY3_1000ms_WNT5A
figure(10)
heatmap_gene10 = heatmap(spotTable,"X","Y","ColorVariable","R3_CY3_1000ms_WNT5A", "ColorMethod","none")
heatmap_gene10.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene10,'R3_CY3_1000ms_WNT5A_heatmap_E158rep2.tif')

% R3_A594_250ms_DDX58
figure(11)
heatmap_gene11 = heatmap(spotTable,"X","Y","ColorVariable","R3_A594_250ms_DDX58", "ColorMethod","none")
heatmap_gene11.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene11,'R3_A594_250ms_DDX58_heatmap_E158rep2.tif')

% R3_CY5_100ms_MITF
figure(12)
heatmap_gene12 = heatmap(spotTable,"X","Y","ColorVariable","R3_CY5_100ms_MITF", "ColorMethod","none")
heatmap_gene12.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene12,'R3_CY5_100ms_MITF_heatmap_E158rep2.tif')

 %% normalize spot counts to housekeeping gene R3_YFP_500ms_UBC **optional**

spotTable.normR1_YFP_500ms_UBC = spotTable.R1_YFP_500ms_UBC./spotTable.R3_YFP_500ms_UBC;
spotTable.normR1_CY3_250ms_NGFR = spotTable.R1_CY3_250ms_NGFR./spotTable.R3_YFP_500ms_UBC;
spotTable.normR1_A594_500ms_MMP1 = spotTable.R1_A594_500ms_MMP1./spotTable.R3_YFP_500ms_UBC;
spotTable.normR1_CY5_250ms_AXL = spotTable.R1_CY5_250ms_AXL./spotTable.R3_YFP_500ms_UBC;
spotTable.normR2_YFP_500ms_UBC = spotTable.R2_YFP_500ms_UBC./spotTable.R3_YFP_500ms_UBC;
spotTable.normR2_CY3_250ms_ITGA3 = spotTable.R2_CY3_250ms_ITGA3./spotTable.R3_YFP_500ms_UBC;
spotTable.normR2_A594_1000ms_FN1 = spotTable.R2_A594_1000ms_FN1./spotTable.R3_YFP_500ms_UBC;
spotTable.normR2_CY5_500ms_EGFR = spotTable.R2_CY5_500ms_EGFR./spotTable.R3_YFP_500ms_UBC;
spotTable.normR3_YFP_500ms_UBC = spotTable.R3_YFP_500ms_UBC./spotTable.R3_YFP_500ms_UBC;
spotTable.normR3_CY3_1000ms_WNT5A = spotTable.R3_CY3_1000ms_WNT5A./spotTable.R3_YFP_500ms_UBC;
spotTable.normR3_A594_250ms_DDX58 = spotTable.R3_A594_250ms_DDX58./spotTable.R3_YFP_500ms_UBC;
spotTable.normR3_CY5_100ms_MITF = spotTable.R3_CY5_100ms_MITF./spotTable.R3_YFP_500ms_UBC;

% check
spotTable(1:50:end,:)

%% make normalized heatmap for each gene

% R1_YFP_500ms_UBC
figure(1)
heatmap_gene1 = heatmap(spotTable,"X","Y","ColorVariable","normR1_YFP_500ms_UBC", "ColorMethod","none")
heatmap_gene1.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene1,'normR1_YFP_500ms_UBC_heatmap_E158rep2.tif')

% R1_CY3_250ms_NGFR
figure(2)
heatmap_gene2 = heatmap(spotTable,"X","Y","ColorVariable","normR1_CY3_250ms_NGFR", "ColorMethod","none")
heatmap_gene2.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene2,'normR1_CY3_250ms_NGFR_heatmap_E158rep2.tif')

% R1_A594_500ms_MMP1
figure(3)
heatmap_gene3 = heatmap(spotTable,"X","Y","ColorVariable","normR1_A594_500ms_MMP1", "ColorMethod","none")
heatmap_gene3.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene3,'normR1_A594_500ms_MMP1_heatmap_E158rep2.tif')

% R1_CY5_250ms_AXL
figure(4)
heatmap_gene4 = heatmap(spotTable,"X","Y","ColorVariable","normR1_CY5_250ms_AXL", "ColorMethod","none")
heatmap_gene4.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene4,'normR1_CY5_250ms_AXL_heatmap_E158rep2.tif')

% R2_YFP_500ms_UBC
figure(5)
heatmap_gene5 = heatmap(spotTable,"X","Y","ColorVariable","normR2_YFP_500ms_UBC", "ColorMethod","none")
heatmap_gene5.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene5,'normR2_YFP_500ms_UBC_heatmap_E158rep2.tif')

% R2_CY3_250ms_ITGA3
figure(6)
heatmap_gene6 = heatmap(spotTable,"X","Y","ColorVariable","normR2_CY3_250ms_ITGA3", "ColorMethod","none")
heatmap_gene6.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene6,'normR2_CY3_250ms_ITGA3_heatmap_E158rep2.tif')

% R2_CY5_500ms_EGFR
figure(7)
heatmap_gene7 = heatmap(spotTable,"X","Y","ColorVariable","normR2_CY5_500ms_EGFR", "ColorMethod","none")
heatmap_gene7.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene7,'normR2_CY5_500ms_EGFR_heatmap_E158rep2.tif')

% R2_A594_1000ms_FN1
figure(8)
heatmap_gene8 = heatmap(spotTable,"X","Y","ColorVariable","normR2_A594_1000ms_FN1", "ColorMethod","none")
heatmap_gene8.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene8,'normR2_A594_1000ms_FN1_heatmap_E158rep2.tif')

% R3_YFP_500ms_UBC
figure(9)
heatmap_gene9 = heatmap(spotTable,"X","Y","ColorVariable","normR3_YFP_500ms_UBC", "ColorMethod","none")
heatmap_gene9.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene9,'normR3_YFP_500ms_UBC_heatmap_E158rep2.tif')

% R3_CY3_1000ms_WNT5A
figure(10)
heatmap_gene10 = heatmap(spotTable,"X","Y","ColorVariable","normR3_CY3_1000ms_WNT5A", "ColorMethod","none")
heatmap_gene10.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene10,'normR3_CY3_1000ms_WNT5A_heatmap_E158rep2.tif')

% R3_A594_250ms_DDX58
figure(11)
heatmap_gene11 = heatmap(spotTable,"X","Y","ColorVariable","normR3_A594_250ms_DDX58", "ColorMethod","none")
heatmap_gene11.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene11,'normR3_A594_250ms_DDX58_heatmap_E158rep2.tif')

% R3_CY5_100ms_MITF
figure(12)
heatmap_gene12 = heatmap(spotTable,"X","Y","ColorVariable","normR3_CY5_100ms_MITF", "ColorMethod","none")
heatmap_gene12.MissingDataColor = [0.8 0.8 0.8]
saveas(heatmap_gene12,'normR3_CY5_100ms_MITF_heatmap_E158rep2.tif')
