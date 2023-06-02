%% Gene count tables of smFISH Data (FateMap) %%
%This script creates gene count tables from the smFISH images of YG102s5.

clear
%
% import combined table for s5
formatSpec = '%d%d%d%d%d%d%d%d%d%d%d%d%C%f';
combinedTable= readtable('FateMap/extractedData/smFISH/YG102s5/Tspots.csv','Format',formatSpec);
% check
head(combinedTable,5);

% define square dimensions in pixels, each square should have 5+ cells
NumberOfPixelsPerSide = 250;

% define square coordinates
combinedTable.SquareCoordX = floor((combinedTable.globalX - min(combinedTable.globalX))/NumberOfPixelsPerSide) + 1;
combinedTable.SquareCoordY = floor((combinedTable.globalY - min(combinedTable.globalY))/NumberOfPixelsPerSide) + 1;
% check
head(combinedTable,5)

%% sum number of spots per gene per square pt1

% set up matrix for xcoord, ycoord, and gene counts
rows = max(combinedTable.SquareCoordX)*max(combinedTable.SquareCoordY); r = 1;
spotNumbers = zeros(rows,5);
for i=1:max(combinedTable.SquareCoordX)
    for j=1:max(combinedTable.SquareCoordY)
        spotNumbers(r,1) = i;
        spotNumbers(r,2) = j;
        r = r+1;
    end
end
% check
spotNumbers(1:100:end,:)

%% Create tables for each gene based on new good spot status

% gene 1: SOX10
% new good status: intensity >= 500, maskID = 0
gene1T = combinedTable(combinedTable.channel == "A594_SOX10" & combinedTable.intensity>499 & combinedTable.maskID == 0,:);
% check
head(gene1T,5)

% gene 2: AXL
% new good status: intensity >= 500, maskID = 0
gene2T = combinedTable(combinedTable.channel == "CY3_AXL" & combinedTable.intensity>499 & combinedTable.maskID == 0,:);
% check
head(gene2T,5)

% gene 3: NGFR
% new good status: intensity >= 800, maskID = 0
gene3T = combinedTable(combinedTable.channel == "CY5_NGFR" & combinedTable.intensity>799 & combinedTable.maskID == 0,:);
% check
head(gene3T,5)

% concatenate good gene spots tables
goodTable = vertcat(gene1T,gene2T,gene3T);
% check
goodTable(1:500:end,:)


%% sum number of spots per gene per square

% fill in gene spot counts in matrix above
r = 1; gene1_sum = 0; gene2_sum = 0; gene3_sum = 0;
for i = 1:max(combinedTable.SquareCoordX)
    for j=1:max(combinedTable.SquareCoordY)
        
        gene1_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel=="A594_SOX10",:);
        gene1_sum = height(gene1_count);
        spotNumbers(r,3) = gene1_sum;
        
    
        gene2_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "CY3_AXL",:);
        gene2_sum = height(gene2_count);
        spotNumbers(r,4) = gene2_sum; 
        
          
        gene3_count = goodTable(goodTable.SquareCoordX==i & goodTable.SquareCoordY==j & goodTable.channel == "CY5_NGFR",:);
        gene3_sum = height(gene3_count);
        spotNumbers(r,5) = gene3_sum;
        
        r = r+1;
    end
end

% check
spotNumbers(1:100:end,:)

%% Add threshold of housekeeping genes to mark square as good

% housekeeping gene: SOX10, anything below threshold will be NaN
housekeeping_threshold = 3;
for i=1:size(spotNumbers,1)
      if spotNumbers(i,3) < housekeeping_threshold
          spotNumbers(i,3) = NaN;
          spotNumbers(i,4) = NaN;
          spotNumbers(i,5) = NaN;
      end
end

% check
spotNumbers(1:100:end,:)

%% Create table and export

% create table
spotTable = array2table(spotNumbers,"VariableNames",["X","Y","SOX10","AXL","NGFR"]);
% check
head(spotTable,5)
% export csv file
writetable(spotTable,'FateMap/extractedData/smFISH/YG102s5/geneCounts_s5.csv')

