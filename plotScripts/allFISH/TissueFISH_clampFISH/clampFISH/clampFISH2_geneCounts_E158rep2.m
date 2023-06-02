%% Gene count tables of clampFISH2.0 Data (FateMap) %%
%This script creates gene count tables from the clampFISH data in E158rep2 tissue subregions.

clear
%
% import table for subregion --> ****change subregion as needed****
formatSpec = '%d%d%d%d%d%d%d%C%f';
regionTable= readtable('spots.csv','Format',formatSpec);
% check
head(regionTable,5);

% define square dimensions in pixels, each square should have 5+ cells
NumberOfPixelsPerSide = 250;

% define square coordinates
regionTable.SquareCoordX = floor((regionTable.x - min(regionTable.x))/NumberOfPixelsPerSide) + 1;
regionTable.SquareCoordY = floor((regionTable.y - min(regionTable.y))/NumberOfPixelsPerSide) + 1;
% check
head(regionTable,5)

%% sum number of spots per gene per square pt1

% set up matrix for xcoord, ycoord, and gene counts
rows = max(regionTable.SquareCoordX)*max(regionTable.SquareCoordY); r = 1;
spotNumbers = zeros(rows,14);
for i=1:max(regionTable.SquareCoordX)
    for j=1:max(regionTable.SquareCoordY)
        spotNumbers(r,1) = i;
        spotNumbers(r,2) = j;
        r = r+1;
    end
end
% check
% spotNumbers(1:10:end,:)

%% Create tables for each gene based on new good spot status

% gene 1: R1_YFP_500ms_UBC
% new good status: intensity >= 150, maskID = 0
gene1T = regionTable(regionTable.channel == "R1_YFP_500ms_UBC" & regionTable.intensity>149 & regionTable.maskID == 0,:);
% check
head(gene1T,5)

% gene 2: R1_CY3_250ms_NGFR
% new good status: intensity >= 50, maskID = 0
gene2T = regionTable(regionTable.channel == "R1_CY3_250ms_NGFR" & regionTable.intensity>49 & regionTable.maskID == 0,:);
% check
head(gene2T,5)

% gene 3: R1_A594_500ms_MMP1
% new good status: intensity >= 60, maskID = 0
gene3T = regionTable(regionTable.channel == "R1_A594_500ms_MMP1" & regionTable.intensity>59 & regionTable.maskID == 0,:);
% check
head(gene3T,5)

% gene 4: R1_CY5_250ms_AXL
% new good status: intensity >= 60, maskID = 0
gene4T = regionTable(regionTable.channel == "R1_CY5_250ms_AXL" & regionTable.intensity>59 & regionTable.maskID == 0,:);
% check
head(gene4T,5)

% gene 5: R2_YFP_500ms_UBC
% new good status: intensity >= 150, maskID = 0
gene5T = regionTable(regionTable.channel == "R2_YFP_500ms_UBC" & regionTable.intensity>149 & regionTable.maskID == 0,:);
% check
head(gene5T,5)

% gene 6: R2_CY3_250ms_ITGA3
% new good status: intensity >= 50, maskID = 0
gene6T = regionTable(regionTable.channel == "R2_CY3_250ms_ITGA3" & regionTable.intensity>49 & regionTable.maskID == 0,:);
% check
head(gene6T,5)

% gene 7: R2_A594_1000ms_FN1
% new good status: intensity >= 100, maskID = 0
gene7T = regionTable(regionTable.channel == "R2_A594_1000ms_FN1" & regionTable.intensity>99 & regionTable.maskID == 0,:);
% check
head(gene7T,5)

% gene 8: R2_CY5_500ms_EGFR
% new good status: intensity >= 80, maskID = 0
gene8T = regionTable(regionTable.channel == "R2_CY5_500ms_EGFR" & regionTable.intensity>79 & regionTable.maskID == 0,:);
% check
head(gene8T,5)

% gene 9: R3_YFP_500ms_UBC
% new good status: intensity >= 150, maskID = 0
gene9T = regionTable(regionTable.channel == "R3_YFP_500ms_UBC" & regionTable.intensity>149 & regionTable.maskID == 0,:);
% check
head(gene9T,5)

% gene 10: R3_CY3_1000ms_WNT5A
% new good status: intensity >= 70, maskID = 0
gene10T = regionTable(regionTable.channel == "R3_CY3_1000ms_WNT5A" & regionTable.intensity>69 & regionTable.maskID == 0,:);
% check
head(gene10T,5)

% gene 11: R3_A594_250ms_DDX58
% new good status: intensity >= 35, maskID = 0
gene11T = regionTable(regionTable.channel == "R3_A594_250ms_DDX58" & regionTable.intensity>34 & regionTable.maskID == 0,:);
% check
head(gene11T,5)

% gene 12: R3_CY5_100ms_MITF
% new good status: intensity >= 35, maskID = 0
gene12T = regionTable(regionTable.channel == "R3_CY5_100ms_MITF" & regionTable.intensity>34 & regionTable.maskID == 0,:);
% check
head(gene12T,5)

% concatenate good gene spots tables
goodTable = vertcat(gene1T,gene2T,gene3T,gene4T,gene5T,gene6T,gene7T,gene8T,gene9T,gene10T,gene11T,gene12T);
% check
goodTable(1:500:end,:)


%% sum number of spots per gene per square

% fill in gene spot counts in matrix above
r = 1; gene1_sum = 0; gene2_sum = 0; gene3_sum = 0; gene4_sum = 0; gene5_sum = 0; gene6_sum = 0; gene7_sum = 0; gene8_sum = 0; gene9_sum = 0; gene10_sum = 0; gene11_sum = 0; gene12_sum = 0;
for i = 1:max(regionTable.SquareCoordX)
    for j=1:max(regionTable.SquareCoordY)
        
        gene1_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel=="R1_YFP_500ms_UBC",:);
        gene1_sum = height(gene1_count);
        spotNumbers(r,3) = gene1_sum;
        
    
        gene2_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "R1_CY3_250ms_NGFR",:);
        gene2_sum = height(gene2_count);
        spotNumbers(r,4) = gene2_sum; 
        
          
        gene3_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "R1_A594_500ms_MMP1",:);
        gene3_sum = height(gene3_count);
        spotNumbers(r,5) = gene3_sum;

        gene4_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel=="R1_CY5_250ms_AXL",:);
        gene4_sum = height(gene4_count);
        spotNumbers(r,6) = gene4_sum;
        
    
        gene5_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "R2_YFP_500ms_UBC",:);
        gene5_sum = height(gene5_count);
        spotNumbers(r,7) = gene5_sum; 
        
          
        gene6_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "R2_CY3_250ms_ITGA3",:);
        gene6_sum = height(gene6_count);
        spotNumbers(r,8) = gene6_sum;

        gene7_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel=="R2_A594_1000ms_FN1",:);
        gene7_sum = height(gene7_count);
        spotNumbers(r,9) = gene7_sum;
        
    
        gene8_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "R2_CY5_500ms_EGFR",:);
        gene8_sum = height(gene8_count);
        spotNumbers(r,10) = gene8_sum; 
        
          
        gene9_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "R3_YFP_500ms_UBC",:);
        gene9_sum = height(gene9_count);
        spotNumbers(r,11) = gene9_sum;

        gene10_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel=="R3_CY3_1000ms_WNT5A",:);
        gene10_sum = height(gene10_count);
        spotNumbers(r,12) = gene10_sum;
        
    
        gene11_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "R3_A594_250ms_DDX58",:);
        gene11_sum = height(gene11_count);
        spotNumbers(r,13) = gene11_sum; 
        
          
        gene12_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "R3_CY5_100ms_MITF",:);
        gene12_sum = height(gene12_count);
        spotNumbers(r,14) = gene12_sum;
        
        r = r+1;
    end
end

% check
spotNumbers(1:10:end,:)

%% Add threshold of housekeeping genes to mark square as good

% housekeeping gene: R3_YFP_500ms_UBC , all spot counts below threshold will be 0 for all genes
housekeeping_threshold = 7;
for i=1:size(spotNumbers,1)
      if spotNumbers(i,11) < housekeeping_threshold
          spotNumbers(i,3) = NaN;
          spotNumbers(i,4) = NaN;
          spotNumbers(i,5) = NaN;
          spotNumbers(i,6) = NaN;
          spotNumbers(i,7) = NaN;
          spotNumbers(i,8) = NaN;
          spotNumbers(i,9) = NaN;
          spotNumbers(i,10) = NaN;
          spotNumbers(i,11) = NaN;
          spotNumbers(i,12) = NaN;
          spotNumbers(i,13) = NaN;
          spotNumbers(i,14) = NaN;
        
      end
end

% check
spotNumbers(1:10:end,:)

%% Create table and export

% create table
spotTable = array2table(spotNumbers,"VariableNames",["X","Y","R1_YFP_500ms_UBC","R1_CY3_250ms_NGFR","R1_A594_500ms_MMP1", "R1_CY5_250ms_AXL", "R2_YFP_500ms_UBC", "R2_CY3_250ms_ITGA3", "R2_A594_1000ms_FN1", "R2_CY5_500ms_EGFR", "R3_YFP_500ms_UBC", "R3_CY3_1000ms_WNT5A", "R3_A594_250ms_DDX58", "R3_CY5_100ms_MITF"]);
% check
head(spotTable,5)
% export csv file ***change as needed****
writetable(spotTable,'/Users/ssa2724/Dropbox (Personal)/FateMap/extractedData/clampFISH2/Tumor_drug/subregion_6_r2c2/subregion6_geneCount2.csv')

