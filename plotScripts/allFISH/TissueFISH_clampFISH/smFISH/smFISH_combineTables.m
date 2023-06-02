%% Combine each subregion's spots.csv tables into a single table
% YG102s1 AND YG102s5
clear


% change these folders to match your local folder structure:
rawDataParentDir='/Users/ssa2724/Dropbox/mouse'; % top-level folder of where the imaging raw data is
extractedDataParentDir='/Users/ssa2724/Dropbox/FateMap';

% inputs used to split up image files into subregions
nrSub=5; % number of subregion rows
ncSub=4; % number of subregion columns

%% smFISH_QC for YG102 s1
% modify these subfolders as appropriate
rawDataDir=fullfile(rawDataParentDir,filesep,'/YG102/s1/largeImageIndividual');
extractedDataDir=fullfile(extractedDataParentDir,filesep,'/extractedData/smFISH/YG102s1');


[subregionFolders,srNums,srRows,srCols]=makeSubregionFolderNames(nrSub,ncSub);

if ~isfolder(extractedDataDir)
    mkdir(extractedDataDir)
end

imgInfo=imfinfo(fullfile(rawDataDir,filesep,subregionFolders{1},filesep,'DAPI.tif')); % use DAPI.tif in first subregion to get subregion pixel rows, columns
imgPixelRows=imgInfo.Height;
imgPixelCols=imgInfo.Width;

Tspots=table();
fprintf('----- getting spots.csv data from %s -----\n',rawDataDir)


for i=1:length(subregionFolders)
    srNum=srNums(i);
    srRow=srRows(i);
    srCol=srCols(i);
    
    subregionFolder=subregionFolders(i);
    fullSubregionFolderPath=fullfile(rawDataDir,filesep,subregionFolder);
    spotsCsvPath=fullfile(fullSubregionFolderPath,'spots.csv');
    
    
    if isfile(spotsCsvPath)
        fprintf('getting on %s\n',subregionFolder)
        TspotsSub=readtable(spotsCsvPath);
        nSpots=height(TspotsSub);
        TspotsSub.srNum=repmat(srNum,nSpots,1);
        TspotsSub.srRow=repmat(srRow,nSpots,1);
        TspotsSub.srCol=repmat(srCol,nSpots,1);

        % Tspots.x is actually pixel row in image
        % Tspots.y is actually pixel column in image
        
        TspotsSub.globalX = TspotsSub.x + imgPixelRows*(srRow-1);
        TspotsSub.globalY = TspotsSub.y + imgPixelCols*(srCol-1);
        TspotsSub=movevars(TspotsSub,{'srNum','srRow','srCol','globalX','globalY'},'Before','spotID'); % move to first few columns
        Tspots=[Tspots;TspotsSub];
    end
end
% output combined Tspots table
writetable(Tspots,fullfile(extractedDataDir,filesep,'Tspots2.csv'))

% optionally plot all the valid points together to check global structure
TspotsShow=Tspots(Tspots.status==1,:); % if status ==1, then spot is above threshold and not masked. This includes all channels though.
TspotsShow=TspotsShow(1:10:end,:); % downsample to help with speed;
cmap=jet(max(srNums));
colorVect=cmap(TspotsShow.srNum,:);
spotSize=2;
figure(1);
scatter(TspotsShow.globalY,TspotsShow.globalX,spotSize,colorVect)
set(gca, 'YDir','reverse')
title('YG102s1')
xlabel('globalY'); xlabel('globalX');

% %% smFISH_QC for YG102 s5
% % modify these subfolders as appropriate
% rawDataDir=fullfile(rawDataParentDir,filesep,'/YG102/s5/largeImageIndividual');
% extractedDataDir=fullfile(extractedDataParentDir,filesep,'/extractedData/smFISH/YG102s5');
% 
% 
% [subregionFolders,srNums,srRows,srCols]=makeSubregionFolderNames(nrSub,ncSub);
% 
% if ~isfolder(extractedDataDir)
%     mkdir(extractedDataDir)
% end
% 
% imgInfo=imfinfo(fullfile(rawDataDir,filesep,subregionFolders{1},filesep,'DAPI.tif')); % use DAPI.tif in first subregion to get subregion pixel rows, columns
% imgPixelRows=imgInfo.Height;
% imgPixelCols=imgInfo.Width;
% 
% Tspots=table();
% fprintf('----- getting spots.csv data from %s -----\n',rawDataDir)
% 
% 
% for i=1:length(subregionFolders)
%     srNum=srNums(i);
%     srRow=srRows(i);
%     srCol=srCols(i);
%     
%     subregionFolder=subregionFolders(i);
%     fullSubregionFolderPath=fullfile(rawDataDir,filesep,subregionFolder);
%     spotsCsvPath=fullfile(fullSubregionFolderPath,'spots.csv');
%     
%     
%     if isfile(spotsCsvPath)
%         fprintf('getting on %s\n',subregionFolder)
%         TspotsSub=readtable(spotsCsvPath);
%         nSpots=height(TspotsSub);
%         TspotsSub.srNum=repmat(srNum,nSpots,1);
%         TspotsSub.srRow=repmat(srRow,nSpots,1);
%         TspotsSub.srCol=repmat(srCol,nSpots,1);
% 
%         % Tspots.x is actually pixel row in image
%         % Tspots.y is actually pixel column in image
%         
%         TspotsSub.globalX = TspotsSub.x + imgPixelRows*(srRow-1);
%         TspotsSub.globalY = TspotsSub.y + imgPixelCols*(srCol-1);
%         TspotsSub=movevars(TspotsSub,{'srNum','srRow','srCol','globalX','globalY'},'Before','spotID'); % move to first few columns
%         Tspots=[Tspots;TspotsSub];
%     end
% end
% % output combined Tspots table
% writetable(Tspots,fullfile(extractedDataDir,filesep,'Tspots2.csv'))
% 
% % optionally plot all the valid points together to check global structure
% TspotsShow=Tspots(Tspots.status==1,:); % if status ==1, then spot is above threshold and not masked. This includes all channels though.
% TspotsShow=TspotsShow(1:10:end,:); % downsample to help with speed;
% cmap=jet(max(srNums));
% colorVect=cmap(TspotsShow.srNum,:);
% spotSize=2;
% figure(2);
% scatter(TspotsShow.globalY,TspotsShow.globalX,spotSize,colorVect)
% set(gca, 'YDir','reverse')
% title('YG102s5')
% xlabel('globalY'); xlabel('globalX');

%% function to generate subregion folder names
function [subregionFolders,srVect,rVect,cVect]=makeSubregionFolderNames(nrSub,ncSub)
% subregion array definition
numSubregions=nrSub*ncSub;
srVect=1:numSubregions;
rVect=reshape(repmat(1:nrSub,ncSub,1),numSubregions,1)';
cVect=repmat(1:ncSub,1,nrSub);
subregionFolders=strcat("subregion",string((srVect)'),"_r",string(rVect)',"_c",string(cVect)');
end

