function guiHandle = launchD2IFGUI(varargin)
    n = inputParser;
    n.addParameter('scanSummary', 'scanSummary.txt', @ischar);
    n.addParameter('imgMasksFile', 'imgMasks.csv', @ischar); 
    n.addParameter('cellMasksFile', 'cellMasks.csv', @ischar); 
    n.addParameter('nucBoundariesFile', 'nucBoundariesIF.csv', @ischar); 
    n.addParameter('cellBoundariesFile', 'cellBoundariesIF.csv', @ischar); 
    n.addParameter('IFquantFile', 'IFquantTable.csv', @ischar); 
    n.addParameter('preStitchedScan', '', @(x) d2utils.checkFile(x));
    n.addParameter('preStitchedScanFilelist','', @(x)validateattributes(x,{'cell'},{'size',[1 nan]}));
    n.addParameter('nucleiMasks', '', @(x) d2utils.checkFile(x));
    n.addParameter('cellMasks', '', @(x) d2utils.checkFile(x));
    n.addParameter('withNuc', 'some', @(x)mustBeMember(x, {'none', 'some', 'all'}));
    n.addParameter('withCyto', true, @islogical);
    n.addParameter('maskResizeFactor', 1, @(x)validateattributes(x,{'numeric'}, {'scalar', '>', 0}));
    n.addParameter('cellPoseTileTable', 'cellPoseTilePositions.csv', @(x) d2utils.checkFile(x));
    
    n.parse(varargin{:});
%----------------------------------------------------------------
%
    if isempty(n.Results.preStitchedScan) && isempty(n.Results.preStitchedScanFilelist)
        if isfile(n.Results.scanSummary)
            scanObj = scanObject('scanSummary', n.Results.scanSummary);
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
            disp('Stitching DAPI channel. This may take a few minutes.')
            scanObj.stitchDAPI();
            disp('Stitching IF channels. This may take a few minutes.')
            scanObj.stitchChannels();
            disp('Saving stitched scans. This may take several minutes.')
            scanObj.saveStitches();
        end
    elseif ~isempty(n.Results.preStitchedScanFilelist)
        fprintf('Loading pre-stitched scans from provided file list (dapi should be first, followed by FISH channels)\n')
        scanObj=scanObject('preStitchedScanFilelist',n.Results.preStitchedScanFilelist);
        scanObj.scanSummaryFile = n.Results.scanSummary;
        if ~isfile(n.Results.scanSummary)
            scanObj.saveScanSummary();
        end
    else
        fprintf('Loading pre-stitched scans.\nThis may take several minutes.\n')
        scanObj = scanObject('scanFile', n.Results.preStitchedScan);
        scanObj.scanSummaryFile = n.Results.scanSummary;
        scanObj.saveScanSummary();
    end
%----------------------------------------------------------------
%
    if isfile(n.Results.imgMasksFile)
        maskObjImg = maskTable(scanObj, n.Results.imgMasksFile);
    else
        fprintf('Unable to detect %s in your current directory. Creating a new masks object\n', n.Results.imgMasksFile)
        maskObjImg = maskTable(scanObj);
        maskObjImg.masksFile = n.Results.imgMasksFile;
    end
%----------------------------------------------------------------
%
    if isfile(n.Results.cellMasksFile)
        maskObjBoundaries = maskTable(scanObj, n.Results.cellMasksFile);
    else
        fprintf('Unable to detect %s in your current directory. Creating a new masks object\n', n.Results.cellMasksFile)
        maskObjBoundaries = maskTable(scanObj);
        maskObjBoundaries.masksFile = n.Results.cellMasksFile;
    end
%----------------------------------------------------------------
%
    if isfile(n.Results.nucBoundariesFile) && isfile(n.Results.cellBoundariesFile) && isfile(n.Results.IFquantFile) %%
        disp('Loading IFboundaries and IFquant objects.')
        IFboundariesObj = d2IF.IFboundaries(scanObj, maskObjBoundaries, n.Results.nucBoundariesFile, n.Results.cellBoundariesFile);
        IFquantObj = d2IF.IFtable(scanObj, maskObjImg, IFboundariesObj, n.Results.IFquantFile);
        IFboundariesObj.addColors();
        IFboundariesObj.addEmptyRows(1000);
        IFquantObj.addEmptyRows(1000);
    else %There may be a more effient way to code the conditions below
        disp('Making new IFboundaries and IFquant objects.')
        IFboundariesObj = d2IF.IFboundaries(scanObj, maskObjBoundaries);
        IFquantObj = d2IF.IFtable(scanObj, maskObjImg, IFboundariesObj);
        IFboundariesObj.nucBoundariesFile = n.Results.nucBoundariesFile;
        IFboundariesObj.cellBoundariesFile = n.Results.cellBoundariesFile;
        IFquantObj.IFquantFile = n.Results.IFquantFile;
        if ~strcmp(n.Results.withNuc, 'none')
            if ~isempty(n.Results.nucleiMasks) %Load pre-stitched cellpose nuclei mask
                fprintf('Loading nuclei mask file:%s.\n', n.Results.nucleiMasks)
%                 IFboundariesObj.loadCellPoseDapi(sprintf('%s_masks.tif', n.Results.cellPoseNuclei), sprintf('%s_outlines.txt',  n.Results.cellPoseNuclei));
                IFboundariesObj.loadCellPoseDapi2(n.Results.nucleiMasks);
                IFboundariesObj.labelMat2nucTable();
            else
                disp('masking dapi')
                if isempty(n.Results.preStitchedScan) && isempty(n.Results.preStitchedScanFilelist)
                    IFboundariesObj.stitchDAPImask();
                else
                     IFboundariesObj.stitchDAPImask2();
                end
                disp('Finding nuclei boundaries.')
                IFboundariesObj.makeNucleiLabelMat();
            end
        end
        
        if n.Results.withCyto
            if ~isempty(n.Results.cellMasks)
                fprintf('Loading cell mask file:%s.\n', n.Results.cellMasks)
                IFboundariesObj.cellPoseCytoFile = n.Results.cellMasks;
%                 IFboundariesObj.loadCellPoseCyto(sprintf('%s_masks.tif', n.Results.cellPoseCyto), sprintf('%s_outlines.txt',  n.Results.cellPoseCyto));
                IFboundariesObj.loadCellPoseCyto2(n.Results.cellMasks);
                IFboundariesObj.labelMat2cytoTable();
                disp('Quantifying nuclei and cytoplasmic IF signal.')
                switch n.Results.withNuc
                    case 'some'
                        IFquantObj.quantAllCytoBoundaries(true);
                    case 'all'
                        IFboundariesObj.assignNucToCyto();
                        IFboundariesObj.addColors();
                        IFquantObj.quantBoundaries();
                    case 'none'
                        IFboundariesObj.makeNucEmpty();
                        IFboundariesObj.addColors();
                        IFquantObj.quantBoundaries();
                end
            elseif ~strcmp(n.Results.withNuc, 'none')
                disp('Quantifying nuclei and cytoplasmic donut.')
                IFquantObj.quantAllLabelMat2();
            else
                disp('If withNuc = none and withCyto = true, then you need to specify cellpose cytoplasmic boundaries')
                return
            end
        else
            if ~strcmp(n.Results.withNuc, 'none')
                %Make empty nuc boundaries
                IFboundariesObj.addColors();
                IFquantObj.quantBoundaries();
            else
                disp('You specified withNuc = none and withCyto = false. Boundaries and quant tables will start out empty')
                IFboundariesObj.addColors();
            end
        end
    end
    IFquantObj.makeCentroidList('meanNuc');
%----------------------------------------------------------------
% 
    disp('Auto-contrasting stitched scans. This may take several minutes.')
    scanObj.contrastDAPIstitch();
    scanObj.contrastStitchedScans([1 99], [0.9 3]);
    disp('Resizing stitched scans')
    scanObj.resizeStitchedScans();
    
    guiHandle = d2IF.d2IFView(scanObj, IFboundariesObj, IFquantObj);
end