function guiHandle = launchD2ThresholdGUI(varargin)
    n = inputParser;
    n.addParameter('scanSummary', 'scanSummary.txt', @ischar);
    n.addParameter('masksFile', 'masks.csv', @ischar); 
    n.addParameter('nucleiFile', 'nuclei.csv', @ischar); 
    n.addParameter('spotsFile', 'spots.csv', @ischar); 
    n.addParameter('preStitchedScan', '', @(x) d2utils.checkFile(x));
    n.addParameter('preStitchedScanFilelist','', @(x)validateattributes(x,{'cell'},{'size',[1 nan]}));
    n.addParameter('channelNames',{}, @(x) isvector(x) && iscellstr(x));
    n.addParameter('channelTypes',{}, @(x) iscell(x) && all(ismember(x,{'dapi','FISH','other'})) && size(x,1)==1);
    n.addParameter('nucleiMasks', '', @(x) d2utils.checkFile(x)); 
    n.addParameter('maskResizeFactor', 1, @(x)validateattributes(x,{'numeric'}, {'vector', '>', 0})); %Consider updating to limit ncol to 1 | 2
    n.addParameter('cellPose', '', @(x) d2utils.checkDir(x));
    n.addParameter('cellPoseTileTable', 'cellPoseTilePositions.csv', @ischar);
    n.addParameter('subtractBackground', false, @islogical);
    n.addParameter('launchGUI',true, @islogical);
    n.addParameter('sigma',[],@(x) validateattributes(x,{'numeric'},{'positive','nonzero'})) % defaults stored in spotsTable.m
    n.addParameter('thresholds',[],@(x) validateattributes(x,{'numeric'},{'size',[1 nan]})) % defaults stored in spotsTable.m
    n.addParameter('aTrousMinThreshFactor', [], @(x)validateattributes(x,{'numeric'}, {'scalar', '>', 0}));
    
    n.parse(varargin{:});
%----------------------------------------------------------------

    if isempty(n.Results.preStitchedScan) && isempty(n.Results.preStitchedScanFilelist)
        if isfile(n.Results.scanSummary)
            scanObj = scanObject('scanSummary', n.Results.scanSummary,'channels', n.Results.channelNames, 'channelTypes',n.Results.channelTypes);
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

        if isfile('stitchedScans.mat')
            fprintf('Loading stitched scans.\nThis may take several minutes.\n')
            scanObj.loadStitches();
        else
            %Need to update this to throw error if no dapi channel
            disp('Stitching DAPI channel. This may take a few minutes.')
            scanObj.stitchDAPI();
            if n.Results.subtractBackground
                disp('Measuring background fluorescence. This may take a few minutes.')
                scanObj.measureBackground;
                disp('Stitching FISH channels. This may take a few minutes.')
                scanObj.stitchChannels2(n.Results.subtractBackground);
            else
                disp('Stitching FISH channels. This may take a few minutes.')
                scanObj.stitchChannels();
            end
            disp('Saving stitched scans. This may take several minutes.')
            scanObj.saveStitches();
        end
        
    elseif ~isempty(n.Results.preStitchedScanFilelist)
        fprintf('Loading pre-stitched scans from provided file list\n')
        scanObj=scanObject('preStitchedScanFilelist',n.Results.preStitchedScanFilelist,'channelTypes',n.Results.channelTypes);
        scanObj.scanSummaryFile = n.Results.scanSummary;
        if ~isfile(n.Results.scanSummary)
            scanObj.saveScanSummary();
        end
    else % then ~isempty(n.Results.preStitchedScan)
        fprintf('Loading pre-stitched scans.\nThis may take several minutes.\n')
        scanObj = scanObject('scanFile', n.Results.preStitchedScan,'channels', n.Results.channelNames,'channelTypes',n.Results.channelTypes);
        %Need to update this to throw error if no dapi channel
        scanObj.scanSummaryFile = n.Results.scanSummary;
        if ~isfile(n.Results.scanSummary)
            scanObj.saveScanSummary();
        end
    end
%----------------------------------------------------------------
%
    if isfile(n.Results.masksFile)
        maskObj = maskTable(scanObj, n.Results.masksFile);
    else
        fprintf('Unable to detect %s in your current directory. Creating a new masks object\n', n.Results.masksFile)
        maskObj = maskTable(scanObj);
        maskObj.masksFile = n.Results.masksFile;
    end
%----------------------------------------------------------------
%
    if isfile(n.Results.nucleiFile)
        disp('Loading nuclei file.')
        nucleiObj = nucleiTable(scanObj, maskObj, n.Results.nucleiFile);
    else
        fprintf('Unable to detect %s in your current directory. Creating a new nuclei object\n', n.Results.nucleiFile)
        nucleiObj = nucleiTable(scanObj, maskObj);
        nucleiObj.nucleiFile = n.Results.nucleiFile;
        if isempty(n.Results.nucleiMasks) && isempty(n.Results.cellPose) %No input nuclei mask
            disp('Finding nuclei. This may take a few minutes.')
            if isempty(n.Results.preStitchedScan) && isempty(n.Results.preStitchedScanFilelist)
                nucleiObj.stitchDAPImask();
            else
                nucleiObj.stitchDAPImask2();
            end
            nucleiObj.findNuclei();
        elseif ~isempty(n.Results.nucleiMasks) %Priority to load pre-stitched masks
            fprintf('Loading nuclei mask file:%s.\nResize factor is %d.\n', n.Results.nucleiMasks, n.Results.maskResizeFactor)
            try %first try loading as image (label matrix)
                nucleiObj.loadLabelMat(n.Results.nucleiMasks, n.Results.maskResizeFactor);
            catch ME
                if (strcmp(ME.identifier,'MATLAB:imagesci:imread:fileFormat'))
                    %try loading as outlines file
                    nucleiObj.loadNucleiOutlines(n.Results.nucleiMasks, n.Results.maskResizeFactor);
                else
                    rethrow(ME)
                end
            end
        elseif ~isempty(n.Results.cellPose) %Load & stitch cellpose masks
            if isfile(n.Results.cellPoseTileTable)
                fprintf('Stitching cellpose masks in directory: %s.\nResize factor is %d.\n', n.Results.cellPose, n.Results.maskResizeFactor)
                nucleiObj.stitchCellPoseMasks(n.Results.cellPoseTileTable, n.Results.cellPose, n.Results.maskResizeFactor);
            else
                fprintf('Unable to find cellpose file table: %s.\nThe file is need for stitching maks in %s\n', n.Results.cellPoseTileTable, n.Results.cellPose)
                disp('If you intended to input a pre-stitched mask, please use the nucleiMasks parameter.')
            end
        end
    end
    nucleiObj.addColors();
    nucleiObj.updateAllMasks();
%----------------------------------------------------------------
    spotsObj = spotTable(scanObj, maskObj, nucleiObj, 'spotsFile',n.Results.spotsFile,'sigma',n.Results.sigma,'thresholds',n.Results.thresholds,'aTrousMinThreshFactor',n.Results.aTrousMinThreshFactor);
    if ~isfile(n.Results.spotsFile)
        disp('Finding spots. This may take several minutes.')
        if spotsObj.minChannelIntensityIsBlockSpecific % Note that if you delete the spotsFile but keep the scanSummary with channel thresholds, it will use these thresholds (and findSpots5) to find spots.  
            spotsObj.findSpots4(); %Using previous findSpots function. 
        else
            spotsObj.findSpots5();
        end
        spotsObj.maskBorderSpots();
        disp('Finished finding spots')
        spotsObj.assignSpotsToNuclei();
    end
    
    if isempty(spotsObj.thresholds) % did not already get thresholds from user input or scanSummary.txt
        spotsObj.defaultThresholds();
    end
    
    spotsObj.updateScanSummary();
    spotsObj.updateAllMasks();
    spotsObj.updateAllSpotStatus();
    spotsObj.makeCentroidList();
    
    %----------------------------------------------------------------
    if ~n.Results.launchGUI % don't launch GUI, just save files
        fprintf('launchGUI==false, so saving mask table, cell table, and spot tables.\nThis may take a minute\n')
        nucleiObj.saveNucleiTable;
        %spotsObj.updateScanSummary;
        spotsObj.saveSpotsTable;
        maskObj.saveMasksTable;
        disp('done')
        
        % fprintf('Saving %s\n',spotsObj.spotsSummaryFile)
        % spotsObj.exportSpotsSummary % could do this by default, but its misleading if thresholds aren't actually good
        % disp('done')
        
        guiHandle=[];
    else % launch GUI
        disp('Auto-contrasting stitched scans. This may take several minutes.')
        scanObj.contrastDAPIstitch();
        scanObj.contrastStitchedScans([1 99], [0.9 3]);
        disp('Resizing stitched scans')
        scanObj.resizeStitchedScans();
        
        guiHandle = d2ThresholdView2(scanObj, maskObj, nucleiObj, spotsObj);
    end
    
end