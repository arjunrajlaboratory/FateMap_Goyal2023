%% Gene count tables of smFISH Data (FateMap) %% 

%This script creates gene count tables from the smFISH data in for resistant colonies (Figure 1). 

clear 

% 

% import spots table *** Note: Adjust file directory to read in table *** 

formatSpec = '%d%d%d%d%d%d%d%C%f'; 

spotsTable= readtable('spots.csv','Format',formatSpec); 

% check 

head(spotsTable,5); 

% define square dimensions in pixels, each square should have 5+ cells 

NumberOfPixelsPerSide = 250; 

% define square coordinates 

spotsTable.SquareCoordX = floor((spotsTable.x - min(spotsTable.x))/NumberOfPixelsPerSide) + 1; 

spotsTable.SquareCoordY = floor((spotsTable.y - min(spotsTable.y))/NumberOfPixelsPerSide) + 1; 

% check 

head(spotsTable,5) 

%% sum number of spots per gene per square pt1 

% set up matrix for xcoord, ycoord, and gene counts 

rows = max(spotsTable.SquareCoordX)*max(spotsTable.SquareCoordY); r = 1; 

spotNumbers = zeros(rows,6); 

for i=1:max(spotsTable.SquareCoordX) 

for j=1:max(spotsTable.SquareCoordY) 

spotNumbers(r,1) = i; 

spotNumbers(r,2) = j; 

r = r+1; 

end 

end 

% check 

spotNumbers(1:100:end,:) 

%% Create tables for each gene based on new good spot status 

% gene 1: ACTA2; channel: a594 

% new good status: intensity >= 159, maskID = 0 

gene1T = spotsTable(spotsTable.channel == "A594" & spotsTable.intensity>158 & spotsTable.maskID == 0,:); 

% check 

head(gene1T,5) 

% gene 2: BGN; channel: CY3 

% new good status: intensity >= 230, maskID = 0 

gene2T = spotsTable(spotsTable.channel == "CY3" & spotsTable.intensity>229 & spotsTable.maskID == 0,:); 

% check 

head(gene2T,5) 

% gene 3: NGFR; Channel: CY5 

% new good status: intensity >= 150, maskID = 0 

gene3T = spotsTable(spotsTable.channel == "CY5" & spotsTable.intensity>149 & spotsTable.maskID == 0,:); 

% check 

head(gene3T,5) 

% background; Channel: YFP 

% new good status: intensity >= 1510, maskID = 0 

gene4T = spotsTable(spotsTable.channel == "YFP" & spotsTable.intensity>1509 & spotsTable.maskID == 0,:); 

% check 

head(gene4T,5) 

% concatenate good gene spots tables 

goodTable = vertcat(gene1T,gene2T,gene3T,gene4T); 

% check 

goodTable(1:500:end,:) 

%% sum number of spots per gene per square 

% fill in gene spot counts in matrix above 

r = 1; sumA594_ACTA2 = 0; sumCY3_BGN = 0; sumCY5_NGFR = 0; sumYFP = 0; 

for i = 1:max(spotsTable.SquareCoordX) 

for j=1:max(spotsTable.SquareCoordY) 

gene1_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel=="A594",:); 

sumA594_ACTA2 = height(gene1_count); 

spotNumbers(r,3) = sumA594_ACTA2; 

gene2_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "CY3",:); 

sumCY3_BGN = height(gene2_count); 

spotNumbers(r,4) = sumCY3_BGN;  

gene3_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "CY5",:); 

sumCY5_NGFR = height(gene3_count); 

spotNumbers(r,5) = sumCY5_NGFR; 

gene4_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "YFP",:); 

sumYFP = height(gene4_count); 

spotNumbers(r,6) = sumYFP; 

r = r+1; 

end 

end 

% check 

spotNumbers(1:100:end,:) 

%% Add threshold of housekeeping genes to mark square as good --> NOT 

% NEEDED 

%  

% % housekeeping gene: UBC, anything below threshold will be NaN 

% housekeeping_threshold = 6; 

% for i=1:size(spotNumbers,1) 

% if spotNumbers(i,5) < housekeeping_threshold 

% spotNumbers(i,3) = NaN; 

% spotNumbers(i,4) = NaN; 

% spotNumbers(i,5) = NaN; 

% end 

% end 

%  

% % check 

% spotNumbers(1:100:end,:) 

%% Create table and export 

% create table 

spotTable = array2table(spotNumbers,"VariableNames",["X","Y","ACTA2","BGN","NGFR", "YFP"]); 

% check 

head(spotTable,5) 

% export csv file 

writetable(spotTable,'geneCounts.csv') 