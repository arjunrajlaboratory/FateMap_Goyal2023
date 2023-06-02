classdef IFboundaries < handle
    
    properties (Access = public)
       scanObj 
       maskObj %Masks for cells and nuclei. Not for IF signal. 
       nucBoundaries %based on cellpose outlines
       nucBoundaries2 %based on labelmat
       cellBoundaries %based on cellpose outlines
       cellBoundaries2 %based on labelmat
       channels
       dapiMask
       dapiLabelMat
       cytoLabelMat
       dapiRP
       cytoRP
       randomColors
       nucBoundariesFile = 'nucBoundariesIF.csv'
       cellBoundariesFile = 'cellBoundariesIF.csv'
       cellPoseNucFile
       cellPoseCytoFile 
       minNucArea = 1000;
       strelRadius = 40;
    end
    
    methods
        function p = IFboundaries(scanObject, maskObj, varargin)
            n = inputParser;
            n.addRequired('scanObject', @(x) isa(x, 'scanObject')) 
            n.addRequired('maskObj', @(x) isa(x, 'maskTable'))
            n.addOptional('nucleiFile', '', @ischar)
            n.addOptional('cellsFile', '', @ischar)
            n.parse(scanObject, maskObj, varargin{:})
            
            p.scanObj = n.Results.scanObject;
            p.maskObj = n.Results.maskObj;
            p.channels = p.scanObj.channels(ismember(p.scanObj.channelTypes,{'FISH'}));
%             p.channels = p.scanObj.stitchedScans.labels;
            p.randomColors = single(d2utils.distinguishable_colors(50));
            if isfile(n.Results.nucleiFile)
                p.nucBoundariesFile = n.Results.nucleiFile;
                tmpNucBoundaries = readtable(p.nucBoundariesFile);
                tmpNucBoundaries = rowfun(@(x) polyshape(x), tmpNucBoundaries, 'InputVariables', {'x', 'y'}, 'SeparateInputs', false, 'GroupingVariables', [{'cellID'}, p.channels], 'OutputVariableNames', 'nucBoundary');
                tmpNucBoundaries(:,'GroupCount') = [];
                tmpNucBoundaries.nucBB = cell2mat(arrayfun(@(x) d2utils.polyshapeBoundingBox(x), tmpNucBoundaries.nucBoundary, 'UniformOutput', false));
                tmpNucBoundaries = convertvars(tmpNucBoundaries, p.channels,'logical');
                p.nucBoundaries2 = movevars(tmpNucBoundaries, {'nucBoundary', 'nucBB'}, 'After', 'cellID');
            else
                cellID = zeros(0,1);
                nucBoundary = repmat(polyshape, 0,1);
                nucBB = single(zeros(0,4));
                p.nucBoundaries2 = table(cellID, nucBoundary, nucBB);
                p.nucBoundaries2(:,p.channels) = num2cell(true(0,numel(p.channels)));
            end
            
            if isfile(n.Results.cellsFile)
                p.cellBoundariesFile = n.Results.cellsFile;
                tmpCellBoundaries = readtable(p.cellBoundariesFile);
                tmpCellBoundaries = rowfun(@(x) polyshape(x), tmpCellBoundaries, 'InputVariables', {'x', 'y'}, 'SeparateInputs', false, 'GroupingVariables', [{'cellID'}, p.channels], 'OutputVariableNames', 'cellBoundary');
                tmpCellBoundaries(:,'GroupCount') = [];
                tmpCellBoundaries.cellBB = cell2mat(arrayfun(@(x) d2utils.polyshapeBoundingBox(x), tmpCellBoundaries.cellBoundary, 'UniformOutput', false));
                tmpCellBoundaries = convertvars(tmpCellBoundaries, p.channels,'logical');
                p.cellBoundaries2 = movevars(tmpCellBoundaries, {'cellBoundary', 'cellBB'}, 'After', 'cellID');
            else
                cellID = zeros(0,1);
                cellBoundary = repmat(polyshape, 0,1);
                cellBB = single(zeros(0,4));
                p.cellBoundaries2 = table(cellID, cellBoundary, cellBB);
                p.cellBoundaries2(:,p.channels) = num2cell(true(0,numel(p.channels)));
            end
        end
        
        function p = addEmptyRows(p, n)
            newNuclei = table('Size', [n, size(p.nucBoundaries2, 2)],...
                'VariableNames', p.nucBoundaries2.Properties.VariableNames,...
                'VariableTypes', varfun(@class, p.nucBoundaries2, 'output','cell'));
            newNuclei{:,p.channels} = false(n, numel(p.channels));
            newNuclei.nucBB = zeros(n, 4);
            if ismember('colors', p.nucBoundaries2.Properties.VariableNames)
                newNuclei.colors = zeros(n, 3);
            end
            p.nucBoundaries2 = [p.nucBoundaries2; newNuclei];
            
            newCells = table('Size', [n, size(p.cellBoundaries2, 2)],...
                'VariableNames', p.cellBoundaries2.Properties.VariableNames,...
                'VariableTypes', varfun(@class, p.cellBoundaries2, 'output','cell'));
            newCells{:,p.channels} = false(n, numel(p.channels));
            newCells.cellBB = zeros(n, 4);
            if ismember('colors', p.cellBoundaries2.Properties.VariableNames)
                newCells.colors = zeros(n, 3);
            end
            p.cellBoundaries2 = [p.cellBoundaries2; newCells];
        end
        
        function p = deleteEmptyRows(p)
%             p.nucBoundaries2(p.nucBoundaries2.cellID == 0, :) = [];
%             p.cellBoundaries2(p.cellBoundaries2.cellID == 0, :) = [];
            %updated function below should remove empty polyshape
            %boundaries
            numNucRegionsArray = [p.nucBoundaries2.nucBoundary.NumRegions];
            p.nucBoundaries2(~numNucRegionsArray, :) = [];
            
            numCellRegionsArray = [p.cellBoundaries2.cellBoundary.NumRegions];
            p.cellBoundaries2(~numCellRegionsArray, :) = [];
        end
        
        function p = makeNucEmpty(p)
            tmpCellIDs = p.cellBoundaries2.cellID;
            tmpBoundaries = repmat({polyshape}, numel(tmpCellIDs), 1);
            tmpBB = zeros(numel(tmpCellIDs), 4);
            tmpStatus = true(numel(tmpCellIDs),numel(p.channels));
            p.nucBoundaries2 = cell2table([num2cell(tmpCellIDs), tmpBoundaries, num2cell(tmpBB, 2), num2cell(tmpStatus)], 'VariableNames', [{'cellID', 'nucBoundary', 'nucBB'}, p.channels]);
        end
        
        function p = makeCytoEmpty(p)
            tmpCellIDs = p.nucBoundaries2.cellID;
            tmpBoundaries = repmat({polyshape}, numel(tmpCellIDs), 1);
            tmpBB = zeros(numel(tmpCellIDs), 4);
            tmpStatus = true(numel(tmpCellIDs),numel(p.channels));
            p.cellBoundaries2 = cell2table([num2cell(tmpCellIDs), tmpBoundaries, num2cell(tmpBB, 2), num2cell(tmpStatus)], 'VariableNames', [{'cellID', 'cellBoundary', 'cellBB'}, p.channels]);
        end
        
        function p = stitchDAPImask(p, varargin)  
            tileTable = p.scanObj.tilesTable;
            tilesTmp = transpose(p.scanObj.scanMatrix);
            tiles = tilesTmp(:);
            height = p.scanObj.tileSize(1);
            width = p.scanObj.tileSize(2);
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'logical');
            channel = find(ismember(p.scanObj.channels,'dapi'));
            reader = bfGetReader(p.scanObj.scanFile);
            iPlane = reader.getIndex(0, channel - 1, 0) + 1;
            
            if nargin == 1
                s = 0.1;
            elseif nargin == 2
                s = varargin{1};
            end
                
            for i = 1:numel(tiles)
                
                reader.setSeries(tiles(i)-1);
                tmpPlane  = bfGetPlane(reader, iPlane);
                
                tileMask = d2utils.makeDAPImask(scale(tmpPlane), 'sensitivity', s);
                
                tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = tileMask;
            end
            reader.close()
            tmpCC = bwconncomp(tmpStitch);
            tmpArea = regionprops(tmpCC, 'area');
            schmutzIdx = [tmpArea.Area] < p.minNucArea;
            schmutz = tmpCC.PixelIdxList(schmutzIdx);
            for i = 1:numel(schmutz)
                tmpStitch(schmutz{i}) = false;
            end
            p.dapiMask = tmpStitch;
        end
        
        function p = stitchDAPImask2(p, varargin) %For pre-stitched scans
            
            n = inputParser;
            n.addParameter('sensitivity', 0.1, @(x)validateattributes(x,{'numeric'}, {'scalar', '>=',0,'<=',1.0}));    
            n.addParameter('blockSize', [1000 1000], @(x)validateattributes(x,{'numeric'}, {'size', [1 2]}));    
            n.parse(varargin{:});
            %Should maybe check that the block size is >[1 1] and < scanDim.
            s = n.Results.sensitivity;
            block = n.Results.blockSize;
            
            function_mask = @(block_struct) imbinarize(scale(block_struct.data),...
                adaptthresh(scale(block_struct.data), s, 'ForegroundPolarity','bright'));
            
            tmpStitch = blockproc(p.scanObj.dapiStitch, block, function_mask, 'BorderSize', [0 0], 'UseParallel', true);
            tmpCC = bwconncomp(tmpStitch);
            tmpArea = regionprops(tmpCC, 'area');
            schmutzIdx = [tmpArea.Area] < p.minNucArea;
            schmutz = tmpCC.PixelIdxList(schmutzIdx);
            for i = 1:numel(schmutz)
                tmpStitch(schmutz{i}) = false;
            end
            p.dapiMask = tmpStitch;
        end
        
        function p = makeNucleiLabelMat(p)
            %First rounds out dapi masks then finds CC and makes
            %label matrix
            tmpBoundaries = bwboundaries(p.dapiMask);
            tmpBoundariesArea = cellfun(@(x) polyarea(x(:,1), x(:,2)), tmpBoundaries);
            tmpBoundaries = tmpBoundaries(tmpBoundariesArea>= p.minNucArea); 
            
            boundaryStitch = zeros(size(p.dapiMask)); %Could make this max of boundaryArray. Or stitchDim
           
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            warning('off', 'MATLAB:polyshape:boundary3Points')
            warning('off', 'MATLAB:polyshape:boolOperationFailed')
            for i = 1:numel(tmpBoundaries) %This could probably be parfor
                [tmpXlim,tmpYlim]  = boundingbox(polyshape(tmpBoundaries{i}));
                tmpMask = poly2mask(tmpBoundaries{i}(:,2)-tmpYlim(1), tmpBoundaries{i}(:,1)-tmpXlim(1), diff(tmpXlim)+1, diff(tmpYlim)+1);
                boundaryStitch(tmpXlim(1):tmpXlim(2), tmpYlim(1):tmpYlim(2)) = imclose(tmpMask, strel('disk', p.strelRadius)); %Consider adjusting operation or strel. 
            end
            dapiCC = bwconncomp(boundaryStitch);
            p.dapiRP = regionprops(dapiCC);
            p.dapiLabelMat = labelmatrix(dapiCC);
            
            tmpBoundaries = bwboundaries(boundaryStitch, 'noholes');
            nucBoundariesArray = cellfun(@(x) polyshape(x), tmpBoundaries, 'UniformOutput', false);
            nucBoundariesArray = cellfun(@(x) union(rmholes(x)), nucBoundariesArray, 'UniformOutput', false); %Remove holes to simplify plotting. May not be necessary.
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), nucBoundariesArray, 'UniformOutput', false);
            status = true(numel(nucBoundariesArray),numel(p.channels));
            p.nucBoundaries2 = cell2table([num2cell((1:numel(nucBoundariesArray))'), nucBoundariesArray, tmpBB, num2cell(status)], 'VariableNames', [{'cellID', 'nucBoundary', 'nucBB'}, p.channels]);
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
            warning('on', 'MATLAB:polyshape:boundary3Points')
            warning('on', 'MATLAB:polyshape:boolOperationFailed')
        end
        
        function p = loadCellPoseDapi(p, labelMatFile, outlineFile)
            tmpLabelMat = imread(labelMatFile);
            if all(size(tmpLabelMat) ==  p.scanObj.stitchDim)
                p.dapiLabelMat =  tmpLabelMat;
                scaleFactor = [1, 1];
            else
                p.dapiLabelMat = imresize(tmpLabelMat, p.scanObj.stitchDim, 'nearest');
                scaleFactor = p.scanObj.stitchDim./size(tmpLabelMat);
            end
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            warning('off', 'MATLAB:polyshape:boundary3Points')
            polyVect = d2utils.parseCellposeOutlines(outlineFile, 'scaleFactor', scaleFactor, 'flip', true);
            polyVect = scale(polyVect, scaleFactor);

            tmpArea = num2cell([polyVect.area]);
            [tmpX, tmpY] = centroid(polyVect);
            tmpCentroid = num2cell([tmpX', tmpY'],2);
            nucBoundariesArray = num2cell(polyVect);
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), nucBoundariesArray, 'UniformOutput', false);
%             tmpArea = cellfun(@(x) polyarea(x(:,1), x(:,2))*scaleFactor,polyArray, 'UniformOutput', false);
%             tmpCentroid = cellfun(@(x) d2utils.poly2centroid(x)*scaleFactor,polyArray, 'UniformOutput', false);
%             tmpBB = cellfun(@(x) d2utils.polygonBoundingBox2(x)*scaleFactor,polyArray, 'UniformOutput', false);
            p.dapiRP = cell2struct([tmpArea', tmpCentroid, tmpBB'], {'Area', 'Centroid', 'BoundingBox'}, 2);
%             nucBoundariesArray = cellfun(@(x) single(x.Vertices), nucBoundariesArray, 'UniformOutput', false);
%             nucBoundariesHeight = cellfun(@(x) size(x, 1), nucBoundariesArray, 'UniformOutput', true);
%             nucIDArray = repelem(1:numel(nucBoundariesHeight), nucBoundariesHeight);
%             p.nucBoundaries = array2table([single(nucIDArray'), vertcat(nucBoundariesArray{:})], 'VariableNames', {'cellID', 'x', 'y'});
            %Try storying boundaries as polyshape.
            status = true(numel(nucBoundariesArray),numel(p.channels));
            p.nucBoundaries2 = cell2table([num2cell((1:numel(nucBoundariesArray))'), nucBoundariesArray', tmpBB', num2cell(status)], 'VariableNames', [{'cellID', 'nucBoundary', 'nucBB'}, p.channels]);
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
            warning('on', 'MATLAB:polyshape:boundary3Points')
        end
        
        function p = loadCellPoseDapi2(p, labelMatFile)
            tmpLabelMat = imread(labelMatFile);
            if all(size(tmpLabelMat) ==  p.scanObj.stitchDim)
                p.dapiLabelMat =  tmpLabelMat;
                %                 scaleFactor = [1, 1];
            else
                p.dapiLabelMat = imresize(tmpLabelMat, p.scanObj.stitchDim, 'nearest');
                %                 scaleFactor = p.scanObj.stitchDim./size(tmpLabelMat);
            end
            p.dapiRP = regionprops(p.dapiLabelMat);
        end
        
        function p = loadCellPoseCyto(p, labelMatFile, outlineFile)
            tmpLabelMat = imread(labelMatFile);
            if all(size(tmpLabelMat) ==  p.scanObj.stitchDim)
                p.cytoLabelMat =  tmpLabelMat;
                scaleFactor = [1, 1];
            else
                p.cytoLabelMat = imresize(tmpLabelMat, p.scanObj.stitchDim, 'nearest');
                scaleFactor = p.scanObj.stitchDim./size(tmpLabelMat);
            end
            
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            warning('off', 'MATLAB:polyshape:boundary3Points')
            polyVect = d2utils.parseCellposeOutlines(outlineFile, 'scaleFactor', scaleFactor, 'flip', true);
            polyVect = scale(polyVect, scaleFactor);

%             tmpArea = num2cell([polyVect.area]);
%             [tmpX, tmpY] = centroid(polyVect);
%             tmpCentroid = num2cell([tmpX', tmpY'],2);
            cellBoundariesArray = num2cell(polyVect);
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), cellBoundariesArray, 'UniformOutput', false);
%             tmpArea = cellfun(@(x) polyarea(x(:,1), x(:,2))*scaleFactor,polyArray, 'UniformOutput', false);
%             tmpCentroid = cellfun(@(x) d2utils.poly2centroid(x)*scaleFactor,polyArray, 'UniformOutput', false);
%             tmpBB = cellfun(@(x) d2utils.polygonBoundingBox2(x)*scaleFactor,polyArray, 'UniformOutput', false);
%             p.dapiRP = cell2struct([tmpArea', tmpCentroid, tmpBB'], {'Area', 'Centroid', 'BoundingBox'}, 2);
            
%             nucBoundariesArray = cellfun(@(x) single(x.Vertices), nucBoundariesArray, 'UniformOutput', false);
%             nucBoundariesHeight = cellfun(@(x) size(x, 1), nucBoundariesArray, 'UniformOutput', true);
%             nucIDArray = repelem(1:numel(nucBoundariesHeight), nucBoundariesHeight);
%             p.nucBoundaries = array2table([single(nucIDArray'), vertcat(nucBoundariesArray{:})], 'VariableNames', {'cellID', 'x', 'y'});
            %Try storying boundaries as polyshape.
            status = true(numel(cellBoundariesArray),numel(p.channels));
            p.cellBoundaries2 = cell2table([num2cell((1:numel(cellBoundariesArray))'), cellBoundariesArray', tmpBB', num2cell(status)], 'VariableNames', [{'cellID', 'cellBoundary', 'cellBB'}, p.channels]);
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
            warning('on', 'MATLAB:polyshape:boundary3Points')
        end
        
         function p = loadCellPoseCyto2(p, labelMatFile)
            tmpLabelMat = imread(labelMatFile);
            if all(size(tmpLabelMat) ==  p.scanObj.stitchDim)
                p.cytoLabelMat =  tmpLabelMat;
%                 scaleFactor = [1, 1];
            else
                p.cytoLabelMat = imresize(tmpLabelMat, p.scanObj.stitchDim, 'nearest');
%                 scaleFactor = p.scanObj.stitchDim./size(tmpLabelMat);
            end
            
            p.cytoRP = regionprops(p.cytoLabelMat);
            
        end
        
        function p = labelMat2nucTable(p)
            nucBoundariesTmp = cell(0, numel(p.dapiRP));
            for i = 1:numel(p.dapiRP) %can make this parfor?
                tmpBB = p.dapiRP(i).BoundingBox;
%                 tmpBB = [fliplr(floor(tmpBB(1:2))), fliplr(tmpBB(3:4))];
                %Create buffer
%                 tmpShift = ceil(tmpBB(4:-1:3)*0.1); previously buffer = 10%, now 10 pixels.
                tmpStart = max([1,1], floor(tmpBB(2:-1:1))-10);
                tmpEnd = min(tmpStart+tmpBB(4:-1:3)+20, p.scanObj.stitchDim); 
                
                tmpRegionMask = p.dapiLabelMat(tmpStart(1):tmpEnd(1), tmpStart(2):tmpEnd(2));
                tmpDapiMask = tmpRegionMask == i;
                tmpNucBoundary = bwboundaries(tmpDapiMask, 'noholes');
                nucBoundariesTmp{i} = cellfun(@(x) x+tmpStart-1, tmpNucBoundary, 'UniformOutput', false);
            end
            %Update nucBoundaries
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            nucBoundariesArray = cellfun(@(x) polyshape(cell2mat(x)), nucBoundariesTmp, 'UniformOutput', false); 
            nucBoundariesArray = cellfun(@(x) union(rmholes(x)), nucBoundariesArray, 'UniformOutput', false); %Remove holes to simplify plotting
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), nucBoundariesArray, 'UniformOutput', false);
            status = true(numel(nucBoundariesArray), numel(p.channels));
            p.nucBoundaries2 = cell2table([num2cell((1:numel(nucBoundariesArray))'), nucBoundariesArray', tmpBB', num2cell(status)], 'VariableNames', [{'cellID', 'nucBoundary', 'nucBB'}, p.channels]);
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
        end
        
        function p = labelMat2cytoTable(p)
            cytoBoundariesTmp = cell(0, size(p.cellBoundaries2, 1));
            for i = 1:size(p.cellBoundaries2, 1) %can make this parfor?
                tmpBB = p.cellBoundaries2.cellBB(i,:);
                %Create buffer
%                 tmpShift = ceil(tmpBB(3:4)*0.1); previously buffer = 10%, now 10 pixels.
                tmpStart = max([1,1], tmpBB(1:2)-10);
                tmpEnd = min(tmpStart+tmpBB(3:4)+20, p.scanObj.stitchDim);
                
                tmpRegionMask = p.cytoLabelMat(tmpStart(1):tmpEnd(1), tmpStart(2):tmpEnd(2));
                tmpCellMask = tmpRegionMask == i;
                tmpCellBoundary = bwboundaries(tmpCellMask, 'noholes');
                cytoBoundariesTmp{i} = cellfun(@(x) x+tmpStart-1, tmpCellBoundary, 'UniformOutput', false);
            end
            %Update cellBoundaries
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            cellBoundariesArray = cellfun(@(x) polyshape(cell2mat(x)), cytoBoundariesTmp, 'UniformOutput', false);
            cellBoundariesArray = cellfun(@(x) union(rmholes(x)), cellBoundariesArray, 'UniformOutput', false); %Remove holes to simplify plotting
            cellBoundariesArray = cellfun(@(x) d2utils.largestPolyRegion(x), cellBoundariesArray, 'UniformOutput', false); %Keep only largest region to simplify plotting. To keep disconnected regions, need to update overlayCells in IFcontroller 
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), cellBoundariesArray, 'UniformOutput', false);
            status = true(numel(cellBoundariesArray), numel(p.channels));
            p.cellBoundaries2 = cell2table([num2cell((1:numel(cellBoundariesArray))'), cellBoundariesArray', tmpBB', num2cell(status)], 'VariableNames', [{'cellID', 'cellBoundary', 'cellBB'}, p.channels]);
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
        end
        
        function p = labelMat2cytoTable2(p)
            cytoBoundariesTmp = cell(0, numel(p.cytoRP));
            for i = 1:numel(p.cytoRP) %can make this parfor?
                tmpBB = p.cytoRP(i).BoundingBox;
                %Create buffer
%                 tmpShift = ceil(tmpBB(4:-1:3)*0.1); previously buffer = 10%, now 10 pixels.
                tmpStart = max([1,1], floor(tmpBB(2:-1:1))-10);
                tmpEnd = min(tmpStart+tmpBB(4:-1:3)+20, p.scanObj.stitchDim);
                
                tmpRegionMask = p.cytoLabelMat(tmpStart(1):tmpEnd(1), tmpStart(2):tmpEnd(2));
                tmpCellMask = tmpRegionMask == i;
                tmpCellBoundary = bwboundaries(tmpCellMask, 'noholes');
                cytoBoundariesTmp{i} = cellfun(@(x) x+tmpStart-1, tmpCellBoundary, 'UniformOutput', false);
            end
            %Update cellBoundaries
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            cellBoundariesArray = cellfun(@(x) polyshape(cell2mat(x)), cytoBoundariesTmp, 'UniformOutput', false);
            cellBoundariesArray = cellfun(@(x) union(rmholes(x)), cellBoundariesArray, 'UniformOutput', false); %Remove holes to simplify plotting
            cellBoundariesArray = cellfun(@(x) d2utils.largestPolyRegion(x), cellBoundariesArray, 'UniformOutput', false); %Keep only largest region to simplify plotting. To keep disconnected regions, need to update overlayCells in IFcontroller 
            tmpBB = cellfun(@(x) d2utils.polyshapeBoundingBox(x), cellBoundariesArray, 'UniformOutput', false);
            status = true(numel(cellBoundariesArray), numel(p.channels));
            p.cellBoundaries2 = cell2table([num2cell((1:numel(cellBoundariesArray))'), cellBoundariesArray', tmpBB', num2cell(status)], 'VariableNames', [{'cellID', 'cellBoundary', 'cellBB'}, p.channels]);
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
        end
        
        function p = assignNucToCyto(p)
            p.deleteEmptyRows; %Delete empty rows
            newNuclei = table('Size', [size(p.cellBoundaries2, 1), size(p.nucBoundaries2, 2)],...
                'VariableNames', p.nucBoundaries2.Properties.VariableNames,...
                'VariableTypes', varfun(@class, p.nucBoundaries2, 'output','cell'));
            newNuclei.cellID = p.cellBoundaries2.cellID;
            extraNuclei = {};
            extraNucleiBB = {};
            for i = 1:size(p.nucBoundaries2, 1)
                tmpCytoInRect = p.getAllCellBoundariesInRect(p.nucBoundaries2.nucBB(i, :));
                overlapIdx = overlaps(p.nucBoundaries2.nucBoundary(i),tmpCytoInRect.cellBoundary);
                overlapN = sum(overlapIdx);
                if overlapN == 0
                    extraNuclei = [extraNuclei, {p.nucBoundaries2.nucBoundary(i)}];
                    extraNucleiBB = [extraNucleiBB, p.nucBoundaries2.nucBB(i,:)];
                elseif overlapN == 1
                    newNucleiIdx = newNuclei.cellID == tmpCytoInRect.cellID(overlapIdx);
                    newNuclei{newNucleiIdx, 'nucBoundary'} = union(newNuclei{newNucleiIdx, 'nucBoundary'}, p.nucBoundaries2.nucBoundary(i));
%                     newNuclei{newNucleiIdx, 'nucBB'} = union(p.nucBoundaries2.nucBB(i));
                elseif overlapN > 1
                    tmpCyto = tmpCytoInRect(overlapIdx, :);
                    for ii = 1:size(tmpCyto, 1)
                        nucInCyto = intersect(tmpCyto.cellBoundary(ii), p.nucBoundaries2.nucBoundary(i));
                        nucInCytoIdx = newNuclei.cellID == tmpCyto.cellID(ii);
                        newNuclei{nucInCytoIdx, 'nucBoundary'} = union(newNuclei{nucInCytoIdx, 'nucBoundary'}, nucInCyto);
                    end
                end
            end
            %Add nucBB
            tmpBB = arrayfun(@(x) d2utils.polyshapeBoundingBox(x), newNuclei.nucBoundary, 'UniformOutput', false);
            newNuclei.nucBB = cell2mat(tmpBB);
            %if there are extraNuclei, add to newNuclei and cellBoundaries2 (with empty polyshape)
            if ~isempty(extraNuclei)
                startCellID = max(newNuclei.cellID) + 1;
                extraNucleiN = numel(extraNuclei);
                status = true(extraNucleiN,numel(p.channels));
                extraNucleiTable = cell2table([num2cell((startCellID:startCellID+extraNucleiN-1)'), extraNuclei', extraNucleiBB', num2cell(status)], 'VariableNames', [{'cellID', 'nucBoundary', 'nucBB'}, p.channels]);
                newNuclei = [newNuclei; extraNucleiTable];
                
                emptyPoly = repmat({polyshape}, extraNucleiN, 1);
                emptyBB = repmat({zeros(1,4)}, extraNucleiN, 1);
                extraCytoTable = cell2table([num2cell((startCellID:startCellID+extraNucleiN-1)'), emptyPoly, emptyBB, num2cell(status)], 'VariableNames', [{'cellID', 'cellBoundary', 'cellBB'}, p.channels]);
                p.cellBoundaries2 = [p.cellBoundaries2; extraCytoTable];
            end
            p.nucBoundaries2 = newNuclei;
        end
        
        function p = removeNucHoles(p)
            numHoles = [p.nucBoundaries2.nucBoundary.NumHoles];
            p.addEmptyRows(sum(numHoles)); %Technically don't need to add rows to cell boundaries.
            tmpNuc = p.nucBoundaries2(numHoles>0, :);
            for i = size(tmpNuc, 1)
                tmpHoles = holes(tmpNuc.nucBoundary(i));
                for ii = 1:numel(p.channels)
                    for iii = numel(tmpHoles)
                        p.maskObj.addMaskLocalCoords(fliplr(tmpHoles(iii).Vertices), p.channels{ii})
                    end
                end
                tmpNuc.nucBoundary(i) = rmholes(tmpNuc.nucBoundary(i));    
            end
        end
        
        function outCellBoundaries = getCellBoundariesInRect(p, channel, rect) %rect specified as [x y nrows ncols]
            
            idx = rectint(p.cellBoundaries2.cellBB,rect)>0 & p.cellBoundaries2{:,channel};
            
            outCellBoundaries = p.cellBoundaries2(idx,:);
        end
        
        function outNucBoundaries = getNucBoundariesInRect(p,channel, rect) %rect specified as [x y nrows ncols]
            
            idx = rectint(p.nucBoundaries2.nucBB,rect)>0 & p.nucBoundaries2{:,channel};
            
            outNucBoundaries = p.nucBoundaries2(idx,:);
        end
        
        function outNucBoundaries = getAllNucBoundariesInRect(p, rect) %rect specified as [x y nrows ncols]
            
            idx = rectint(p.nucBoundaries2.nucBB,rect)>0;
            
            outNucBoundaries = p.nucBoundaries2(idx,:);
        end
        
        function outCellBoundaries = getAllCellBoundariesInRect(p, rect) %rect specified as [x y nrows ncols]
            
            idx = rectint(p.cellBoundaries2.cellBB,rect)>0;
            
            outCellBoundaries = p.cellBoundaries2(idx,:);
        end
        
        function p = addColors(p)
            %Note: assigns colors to all cells and nuclei, even if status = false. 
            colorArray = [repmat(p.randomColors, floor(size(p.nucBoundaries2, 1)/size(p.randomColors, 1)), 1);...
                p.randomColors(1:mod(size(p.nucBoundaries2, 1), size(p.randomColors, 1)), :)];
           
            p.nucBoundaries2.colors = colorArray;
            [~, cellIdx] = ismember(p.nucBoundaries2.cellID, p.cellBoundaries2.cellID); %in case there is variable # or order of cellID (there shouldn't be). 
            p.cellBoundaries2.colors = single(zeros(numel(p.cellBoundaries2.cellID), 3));
            p.cellBoundaries2.colors(cellIdx,:) = colorArray;
        end
        
        function p = addNewCell(p, newCellID, nucPoly, cellPoly)
            nucBB = d2utils.polyshapeBoundingBox(nucPoly);
            cellBB = d2utils.polyshapeBoundingBox(cellPoly);
            idx = find(p.nucBoundaries2.cellID == 0, 1, 'first');
            randomColor = single(p.randomColors(randi(50),:));
            p.nucBoundaries2(idx,:) = [{newCellID, nucPoly, nucBB}, num2cell(true(1, numel(p.channels))), {randomColor}];
            p.cellBoundaries2(idx,:) = [{newCellID, cellPoly, cellBB}, num2cell(true(1, numel(p.channels))), {randomColor}];
        end
        
        function p = updateNucPoly(p, cellID, nucPoly)
            nucBB = d2utils.polyshapeBoundingBox(nucPoly);
            nucBoundaryIdx = p.nucBoundaries2.cellID == cellID; 
            p.nucBoundaries2(nucBoundaryIdx, {'nucBoundary', 'nucBB'}) = {nucPoly, nucBB};
            p.nucBoundaries2{nucBoundaryIdx, p.channels} = true(1, numel(p.channels));
        end
        
        function p = updateCellPoly(p, cellID, cellPoly)
            cellBB = d2utils.polyshapeBoundingBox(cellPoly);
            cellBoundaryIdx = p.cellBoundaries2.cellID == cellID; 
            p.cellBoundaries2(cellBoundaryIdx, {'cellBoundary', 'cellBB'}) = {cellPoly, cellBB};
            p.cellBoundaries2{cellBoundaryIdx, p.channels} = true(1, numel(p.channels));
        end
        
        function cellIDToMask = addCellMask(p, channel, maskPosition)
            tmpBB = d2utils.polygonBoundingBox(fliplr(maskPosition));
            tmpMaskPoly = polyshape(fliplr(maskPosition));
            
            nucTableTmp = p.getNucBoundariesInRect(channel, tmpBB);
            cellTableTmp = p.getCellBoundariesInRect(channel, tmpBB);
            
            nucIDs = nucTableTmp.cellID(overlaps(tmpMaskPoly, nucTableTmp.nucBoundary));
            cellIDs = cellTableTmp.cellID(overlaps(tmpMaskPoly, cellTableTmp.cellBoundary));
            cellIDToMask = union(nucIDs, cellIDs);
            
            nucIdx = ismember(p.nucBoundaries2.cellID, cellIDToMask);
            p.nucBoundaries2{nucIdx, channel} = false;
            
            cellIdx = ismember(p.cellBoundaries2.cellID, cellIDToMask);
            p.cellBoundaries2{cellIdx, channel} = false;
        end
        
        function p = saveBoundaries(p)
            if isempty(p.nucBoundaries2) || ~any(p.nucBoundaries2{:,p.channels}, 'all')
                fprintf("There are no valid nuclei boundaries to save")
            else
                outTable = p.nucBoundaries2(~p.nucBoundaries2.cellID == 0,:); %Remove empty rows but keep masked cells (status = false). 
                outTable = rowfun(@(x) x.Vertices, outTable, 'InputVariables', 'nucBoundary', 'GroupingVariables', [{'cellID'}, p.channels], 'OutputVariableNames', 'Vertices');
                outTable(:,'GroupCount') = [];
                outTable = splitvars(outTable, 'Vertices', 'NewVariableNames', {'x', 'y'});
                writetable(outTable, p.nucBoundariesFile);
            end
            
            if isempty(p.cellBoundaries2) || ~any(p.cellBoundaries2{:,p.channels}, 'all')
                fprintf("There are no valid cell boundaries to save")
            else
                outTable = p.cellBoundaries2(~p.cellBoundaries2.cellID == 0,:); %Remove empty rows but keep masked cells (status = false). 
                outTable = rowfun(@(x) x.Vertices, outTable, 'InputVariables', 'cellBoundary', 'GroupingVariables', [{'cellID'}, p.channels], 'OutputVariableNames', 'Vertices');
                outTable(:,'GroupCount') = [];
                outTable = splitvars(outTable, 'Vertices', 'NewVariableNames', {'x', 'y'});
                writetable(outTable, p.cellBoundariesFile);
            end
        end

%         function p = maskCellsInChannel(p, channel)
%             maxCellMask = max(p.maskObj.masksBB{p.maskObj.masksBB.dapi,'maskID'});
%             maskBB = p.maskObj.masksBB{p.maskObj.masksBB.maskID == maxCellMask,'BB'}; %Only query nuclei within mask bouding box
%             tmpMaskPoly = polyshape(p.maskObj.masks{p.maskObj.masks.maskID == maxCellMask,{'x', 'y'}});
%             nucTableTmp = p.getNucBoundariesInRect(channel, maskBB);
%             cellTableTmp = p.getCellBoundariesInRect(channel, maskBB);
%             
%             nucIdx = p.getNucleiInRectIdx(maskBB);
%             polyIdx = inpolygon(p.nuclei.x(nucIdx), p.nuclei.y(nucIdx),...
%                 p.maskObj.masks{p.maskObj.masks.maskID == maxCellMask, 'x'}, p.maskObj.masks{p.maskObj.masks.maskID == maxCellMask, 'y'});
%             
%             nucIdx(nucIdx) = polyIdx; %Only nuclei in polygon remain true
%             p.nuclei.maskID(nucIdx) = maxCellMask;
%             p.nuclei.status(nucIdx) = false;
%             
%         end
%         
%         function p = addMask(p, channel, rect, maskCoords)
%             tmpMaskPoly = polyshape(maskCoords);
%             nucTableTmp = p.getNucBoundariesInRect(channel, rect);
%             cellTableTmp = p.getCellBoundariesInRect(channel, rect);
%             
%             nucIDs = nucTableTmp.cellID(overlaps(tmpMaskPoly, nucTableTmp.nucBoundary));
%             cellIDs = cellTableTmp.cellID(overlaps(tmpMaskPoly, cellTableTmp.cellBoundary));
%             cellIDToMask = union(nucIDs, cellIDs);
%             mask
%                 tmpCellIDs = union(nucTableTmp.cellID, cellTableTmp.cellID);
%         end
    end
end