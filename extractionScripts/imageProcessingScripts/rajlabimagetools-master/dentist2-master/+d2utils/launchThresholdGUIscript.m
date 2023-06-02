scanSummaryFile = 'scanSummary.txt';
masksFile = 'masks.csv';
nucleiFile = 'nuclei.csv';
spotsFile = 'spots.csv';
stitchedScansFile = 'stitchedScans.mat';
%----------------------------------------------------------------
%
if isfile(scanSummaryFile)
    scanObj = scanObject('scanSummary', scanSummaryFile);
else
    fprintf('Unable to detect %s in your current directory.\n. Make sure to run the d2StitchingGUI before launching the d2ThresholdGUI.\n You may also want to check your path and %s and try again. ', n.Results.scanSummary, n.Results.scanSummary)
    return
end
    
%Check if stitches have been saved. If not, stitch and save to default
%files. 
if isempty(scanObj.tilesTable)
    disp('The scan object does not contain a tiles table. Creating a new tiles table.')
    scanObj.loadTiles();
    scanObj.savetilesTable();
end

if isfile(stitchedScansFile)
    fprintf('Loading stitched scans.\nThis may take several minutes.\n')
    scanObj.loadStitches();
else
    disp('Stitching DAPI channel. This may take a few minutes.')
    scanObj.stitchDAPI();
    disp('Stitching FISH channels. This may take a few minutes.')
    scanObj.stitchChannels();
    disp('Saving stitched scans. This may take several minutes.')
    scanObj.saveStitches();
end
   
%----------------------------------------------------------------
%   
if isfile(masksFile)
    maskObj = maskTable(scanObj, masksFile);
else
    fprintf('Unable to detect %s in your current directory. Creating a new masks object\n', masksFile)
    maskObj = maskTable(scanObj);
end
%----------------------------------------------------------------
%   
if isfile(nucleiFile)
    disp('Loading nuclei file.')
    nucleiObj = nucleiTable(scanObj, maskObj, nucleiFile);
    nucleiObj.addColors();
else
    fprintf('Unable to detect %s in your current directory. Creating a new nuclei object\n', nucleiFile)
    nucleiObj = nucleiTable(scanObj, maskObj);
    disp('Finding nuclei. This may take a few minutes.')
    nucleiObj.stitchDAPImask();
    nucleiObj.findNuclei();
    nucleiObj.addColors();
    disp('Saving nuclei file.')
    nucleiObj.saveNucleiTable();
end
nucleiObj.updateAllMasks();    
%----------------------------------------------------------------
% 
if isfile(spotsFile)
    spotsObj = spotTable(scanObj, maskObj, nucleiObj, scanSummaryFile, spotsFile);
else
    fprintf('Unable to find %s in your current directory. Creating a new spots object\n', spotsFile)
    spotsObj = spotTable(scanObj, maskObj, nucleiObj, scanSummaryFile);
    disp('Finding spots. This may take several minutes.')
    spotsObj.findSpots3();
    spotsObj.maskBorderSpots();
    disp('Finished finding spots')
    spotsObj.assignSpotsToNuclei();
end

if isempty(spotsObj.thresholds)
    spotsObj.defaultThresholds();
end

spotsObj.updateScanSummary();
spotsObj.updateAllMasks();
spotsObj.updateAllSpotStatus();
spotsObj.makeCentroidList();

disp('Auto-contrasting stitched scans. This may take several minutes.')
scanObj.contrastDAPIstitch();
scanObj.contrastStitchedScans([1 99], [0.9 3]);
disp('Resizing stitched scans')
scanObj.resizeStitchedScans();



guiHandle2 = d2ThresholdView2(scanObj, maskObj, nucleiObj, spotsObj);