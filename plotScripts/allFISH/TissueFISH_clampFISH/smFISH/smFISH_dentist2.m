%% FateMap smFISH analysis
%   Ran this analysis using MATLAB 2022a on an M1 Mac
% 
% Dependencies
%   
%   dentist2
%       https://github.com/arjunrajlaboratory/dentist2
%       version used for this analysis: b1e4b14 on Dec 14, 2021
% 
%   rajlabimagetools
%       https://github.com/arjunrajlaboratory/rajlabimagetools
%       version used for this analysis: 89e9d9f on Sep 17, 2020
% 
%   MATLAB Toolboxes needed:
% 
%     {'Deep Learning Toolbox'                  }
%     {'Image Processing Toolbox'               }
%     {'Statistics and Machine Learning Toolbox'}
%     {'Curve Fitting Toolbox'                  }
%     {'Parallel Computing Toolbox'             }
%     {'MATLAB Parallel Server'                 }
%     {'Polyspace Bug Finder'                   }
%
%
% Datasets include:
% 
% Channel:  CY3     A594        CY5         CY7
% Dye:      Cy3     Alexa594    Atto647N    Atto680
% YG102 s1  BGN     ACTA2       UBC
% YG102 s5  AXL     SOX10       NGFR
% YG106 s4  ACTA2   AXL         UBC         GAPDH   


%% Section 1: general inputs
% change these 3 folders to match your local folder structure:
clear
addpath('/Users/iandardani/Code/dentist2') % dentist2
addpath('/Users/iandardani/Code/rajlabimagetools') % if don't have this, wil get this error: Undefined function 'scale' for input arguments of type 'uint16'. Note that it includes a function 'scale.m' which is also the name of functions in various MATLAB toolboxes
parentDir='/Users/iandardani/Dropbox (RajLab)/mouse'; % top-level folder of where the imaging raw data is

% inputs used to split up image files into subregions
nrSub=5; % number of subregion rows
ncSub=4; % number of subregion columns

%% Section 2: subdivide lage stitched images into smaller 'subregion' images. Dentist2/some computers cannot handle such large images
% % this part can be skipped if the images are already in th 'subregion' folders
% 
% YG102 s1
% Channel   CY3     A594        CY5         CY7
% Dye       Cy3     Alexa594    Atto647N    Atto680
% gene      BGN     ACTA2       UBC
rawDataSubDir='YG102/s1/largeImageIndividual';
subregionFolders=makeSubregionFolderNames(nrSub,ncSub);
preStitchedScanFilelist={...
    'DAPI.tif',...
    'CY3_BGN.tif',...
    'A594_ACTA2.tif',...
    'CY5_UBC.tif',...
    'Brightfield.tif'};
rawDataDir=fullfile(parentDir,filesep,rawDataSubDir);
splitImgIntoSubregions(rawDataDir,preStitchedScanFilelist,nrSub,ncSub)

% YG102 s5
% Channel   CY3     A594        CY5         CY7
% Dye       Cy3     Alexa594    Atto647N    Atto680
% gene      AXL     SOX10       NGFR
rawDataSubDir='YG102/s5/largeImageIndividual';
subregionFolders=makeSubregionFolderNames(nrSub,ncSub);
preStitchedScanFilelist={...
    'DAPI.tif',...
    'A594_SOX10.tif',...
    'CY3_AXL.tif',...
    'CY5_NGFR.tif',...
    'Brightfield.tif'};
rawDataDir=fullfile(parentDir,filesep,rawDataSubDir);
splitImgIntoSubregions(rawDataDir,preStitchedScanFilelist,nrSub,ncSub)

%% Section 3: Process 'subregion' images with dentist2 to find spots
sigma=1; % theoretically, the point spread function at 60X/1.4NA with 2x2 camera binning should have sigma of 0.4. Go larger than this. dentist2's default is 2.
aTrousMinThreshFactor=1.5; % only output spots into spots.csv that are 1.5-fold lower than the threshold (if one is provided), otherwise 1.5-fold lower than the autothreshold for a given block). Eg. if threshold provided is 45, then every spots 30 or greater will be in spots.csv, although all spots <45 will have valid=false
launchGUI=false; % false to run a continuous loop to process spots for all subregions in a batch. Then, afterwards, you can turn to true to QC check these. The second time you run launchD2ThresholdGUI there will be a spots.csv table in the folder, and it will take these instead of finding spots again

% YG102 s1
% Channel   CY3     A594        CY5         CY7
% Dye       Cy3     Alexa594    Atto647N    Atto680
% gene      BGN     ACTA2       UBC
rawDataSubDir='YG102/s1/largeImageIndividual';
subregionFolders=makeSubregionFolderNames(nrSub,ncSub);
preStitchedScanFilelist={...
    'DAPI.tif',...
    'CY3_BGN.tif',...
    'A594_ACTA2.tif',...
    'CY5_UBC.tif',...
    'Brightfield.tif'};
channelTypes={'dapi','FISH','FISH','FISH','other'};
thresholds=   [       500     500   800];
for iSub=1:length(subregionFolders)
    subregionFolder=subregionFolders(iSub);
    fprintf('-------------- working on %s --------------\n',subregionFolder)
    cd(fullfile(parentDir,filesep,rawDataSubDir,filesep,subregionFolder))
    h=launchD2ThresholdGUI('preStitchedScanFilelist',preStitchedScanFilelist,'launchGUI',launchGUI,'sigma',sigma,'channelTypes',channelTypes,'thresholds',thresholds,'aTrousMinThreshFactor',aTrousMinThreshFactor);
    pause(2)
    %clear h
end

%% YG102 s5  
% Channel   CY3     A594        CY5         CY7
% Dye       Cy3     Alexa594    Atto647N    Atto680
% gene      AXL     SOX10       NGFR
rawDataSubDir='YG102/s5/largeImageIndividual';
subregionFolders=makeSubregionFolderNames(nrSub,ncSub);
% some subregions have no cells, exclude those from the list:
subregionFolders=subregionFolders(~ismember(subregionFolders,["subregion1_r1_c1","subregion2_r1_c2","subregion5_r2_c1"]));
preStitchedScanFilelist={...
    'DAPI.tif',...
    'CY3_AXL.tif',...
    'A594_SOX10.tif',...
    'CY5_NGFR.tif',...
    'Brightfield.tif'};
channelTypes={'dapi','FISH','FISH','FISH','other'};
thresholds=   [       700     400   500]; % chosen 5-Sep-2022 2:45pm based on iSub=4
for iSub=1:length(subregionFolders)
    subregionFolder=subregionFolders(iSub);
    fprintf('-------------- working on %s --------------\n',subregionFolder)
    cd(fullfile(parentDir,filesep,rawDataSubDir,filesep,subregionFolder))
    h=launchD2ThresholdGUI('preStitchedScanFilelist',preStitchedScanFilelist,'launchGUI',launchGUI,'sigma',sigma,'channelTypes',channelTypes,'thresholds',thresholds,'aTrousMinThreshFactor',aTrousMinThreshFactor);
    pause(2)
    %clear h
end


%% functions
%  function to output the subregion images into their subfolders
function splitImgIntoSubregions(imgDir,fileList,nrSub,ncSub)
numSubregions=nrSub*ncSub;
origDir=pwd;
cd(imgDir)

subregionFolders=makeSubregionFolderNames(nrSub,ncSub);
srVect=1:(nrSub*ncSub);
% remove current subregion folders
for iSub=srVect
    subregionFolder=subregionFolders(iSub);
    if isfolder(subregionFolder)
        rmdir(subregionFolder,"s")
    end
end
% check if there are other subregion folders
temp=dir;
filesFolders=[{temp.name}];
if any(all([startsWith(filesFolders,'subregion','IgnoreCase',true)',[temp.isdir]'],2))
    error("clear folders starting with 'subregion' in %s",imgDir)
end
for iSub=srVect
    subregionFolder=subregionFolders(iSub);
    mkdir(subregionFolder)
end
% determine row and column start/end
finfo=imfinfo(fileList{1});
ncRaw=finfo.Width;
nrRaw=finfo.Height;
%imgRaw=imread(fileList{1});
%[nrRaw,ncRaw]=size(imgRaw);
nrImg=floor(nrRaw/nrSub); % number of pixel rows in subregion image
ncImg=floor(ncRaw/ncSub); % number of pixel columns in subregion image
rStarts1=1:nrImg:(nrImg*nrSub);
cStarts1=1:ncImg:(ncImg*ncSub);
rStarts=reshape(repmat(rStarts1,ncSub,1),numSubregions,1)';
cStarts=repmat(cStarts1,1,nrSub);
rEnds=rStarts+nrImg-1;
cEnds=cStarts+ncImg-1;

% crop each subregion & write to image
for iFile=1:length(fileList)
    fileName=fileList{iFile};
    imgRaw=imread(fileName); % will read first plane if multi-plane
    assert(isequal(size(imgRaw),[nrRaw,ncRaw]))

    for iSub=srVect
        subregionFolder=subregionFolders(iSub);
        outName=fullfile(subregionFolder,filesep,fileName);
        imgSub=imgRaw(rStarts(iSub):rEnds(iSub),cStarts(iSub):cEnds(iSub));
        imwrite(imgSub,outName)
    end
end


cd(origDir)
end

% function to generate subregion folder names
function [subregionFolders,srVect]=makeSubregionFolderNames(nrSub,ncSub)
% subregion array definition
numSubregions=nrSub*ncSub;
srVect=1:numSubregions;
rVect=reshape(repmat(1:nrSub,ncSub,1),numSubregions,1)';
cVect=repmat(1:ncSub,1,nrSub);
subregionFolders=strcat("subregion",string((srVect)'),"_r",string(rVect)',"_c",string(cVect)');
end
