% smFISH_QC for YG102 s5
clear

% change these 2 folders to match your local folder structure:
addpath('/Users/ssa2724/Documents/GitHub/dentist2/') % dentist2
addpath('/Users/ssa2724/Documents/GitHub/rajlabimagetools/') % if don't have this, wil get this error: Undefined function 'scale' for input arguments of type 'uint16'. Note that it includes a function 'scale.m' which is also the name of functions in various MATLAB toolboxes



sigma=1; % theoretically, the point spread function at 60X/1.4NA with 2x2 camera binning should have sigma of 0.4. Go larger than this. dentist2's default is 2.
aTrousMinThreshFactor=1.5; % only output spots into spots.csv that are 1.5-fold lower than the threshold (if one is provided), otherwise 1.5-fold lower than the autothreshold for a given block). Eg. if threshold provided is 45, then every spots 30 or greater will be in spots.csv, although all spots <45 will have valid=false
launchGUI=true; % false to run a continuous loop to process spots for all subregions in a batch. Then, afterwards, you can turn to true to QC check these. The second time you run launchD2ThresholdGUI there will be a spots.csv table in the folder, and it will take these instead of finding spots again

% YG102 s5
% Channel   CY3     A594        CY5         CY7
% Dye       Cy3     Alexa594    Atto647N    Atto680
% gene      AXL     SOX10       NGFR
preStitchedScanFilelist={...
    'DAPI.tif',...
    'A594_SOX10.tif',...
    'CY3_AXL.tif',...
    'CY5_NGFR.tif',...
    'Brightfield.tif'};
channelTypes={'dapi','FISH','FISH','FISH','other'};
thresholds=   [       500     500   800];
h=launchD2ThresholdGUI('preStitchedScanFilelist',preStitchedScanFilelist,'launchGUI',launchGUI,'sigma',sigma,'channelTypes',channelTypes,'thresholds',thresholds,'aTrousMinThreshFactor',aTrousMinThreshFactor);

