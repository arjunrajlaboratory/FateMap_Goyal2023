classdef nucleiTable < handle
    
    properties (Access = public)
        
        scanObj
        maskObj
        
        nuclei
        minNucleusSize = 1000; %Update set method to call findNuclei and updateAllMasks
        nucMaskSens = 0.1;
        dapiMask
        nucleiFile = 'nuclei.csv'
        nucleiChanged = false %kinda serving as an event
    end
    
    events
        %Might be worth making events for when cell table changes
        %(add/mask/remove cells), and then add listeners in mainAxesCntrlr.
    end
    
    methods
        
        function p = nucleiTable(scanObject, maskObj, varargin)
            p.scanObj = scanObject;
            p.maskObj = maskObj;
            if nargin == 2
                fprintf('New Table\n');
                p.nuclei = table('size', [0,7],... %Possibly unnecessary since new table created with p.findNuclei. 
                    'VariableNames', {'nucID', 'x', 'y', 'status', 'maskID', 'nucleusArea', 'colors'},...
                    'VariableTypes', [repmat({'single'}, 1, 3), {'logical'}, repmat({'single'}, 1, 3)]);
            elseif nargin == 3
                fprintf('Loading Table\n');
                p.nucleiFile = varargin{1};
                opts = detectImportOptions(varargin{1});
                opts = setvartype(opts, 'single');
                p.nuclei = readtable(varargin{1},opts);
                p.nuclei.status = logical(p.nuclei.status); %For some reason, when I set 'status' to 'logical' they all go to false. So doing this instead
            end
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
                p.nucMaskSens = 0.1;
            elseif nargin == 2
                p.nucMaskSens = varargin{1};
            end
                
            for i = 1:numel(tiles)
                
                reader.setSeries(tiles(i)-1);
                tmpPlane  = bfGetPlane(reader, iPlane);
                
                tileMask = d2utils.makeDAPImask(scale(tmpPlane), 'sensitivity', p.nucMaskSens);
                
                tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = tileMask;
            end
            reader.close()
            
            p.dapiMask = tmpStitch;
        end
        
        function p = stitchDAPImask2(p, varargin) %For pre-stitched scans
            
            n = inputParser;
            n.addParameter('sensitivity', 0.1, @(x)validateattributes(x,{'numeric'}, {'scalar', '>=',0,'<=',1.0}));    
            n.addParameter('blockSize', [1500 1500], @(x)validateattributes(x,{'numeric'}, {'size', [1 2]}));    
            n.parse(varargin{:});
            %Should maybe check that the block size is >[1 1] and < scanDim.
            p.nucMaskSens = n.Results.sensitivity;
            block = n.Results.blockSize;
            
            function_mask = @(block_struct) imbinarize(scale(block_struct.data),...
                adaptthresh(scale(block_struct.data), p.nucMaskSens, 'ForegroundPolarity','bright'));
            
            p.dapiMask = blockproc(p.scanObj.dapiStitch, block, function_mask, 'BorderSize', [0 0], 'UseParallel', true);
        end
        
        function p = findNuclei(p)
            
            CC = bwconncomp(p.dapiMask); %May be able to delete this and run regionProps on dapiMask
            rp = regionprops(CC);
            area = [rp.Area];
            idx = area >= p.minNucleusSize;
            rp = rp(idx);
            centroids = [rp.Centroid];
            centroids = round(reshape(centroids,2,[])');
            centroids = single(centroids);
            area = single([rp.Area]);
            
            numCentroids = size(centroids, 1);
            status = true(numCentroids, 1);
            maskID = single(zeros(numCentroids, 1));
            colors = single(zeros(numCentroids,3));

            p.nuclei = table((1:numCentroids)',centroids(:,2),centroids(:,1), status, maskID, area', colors,...
                'VariableNames', {'nucID', 'x', 'y', 'status', 'maskID', 'nucleusArea', 'colors'});
            p.addEmptyRows(1000);
        end
        
        function p = loadNucleiOutlines(p, inFileName, scaleFactor)
            warning('off', 'MATLAB:polyshape:repairedBySimplify') 
            polyArray = d2utils.parseCellposeOutlines(inFileName);
            warning('on', 'MATLAB:polyshape:repairedBySimplify')
            if isempty(polyArray)
                error('Unable to parse masks from %s. Please check file format.', inFileName)
            end
            polyArray = scale(polyArray, scaleFactor);
            [polyX, polyY] = centroid(polyArray);
            polyArea = single(area(polyArray));
            
            status = true(numel(polyX), 1);
            maskID = single(zeros(numel(polyX), 1));
            colors = single(zeros(numel(polyX),3));

            p.nuclei = table((1:numel(polyX))',single(polyY)',single(polyX)', status, maskID, polyArea', colors,...
                'VariableNames', {'nucID', 'x', 'y', 'status', 'maskID', 'nucleusArea', 'colors'});
            p.addEmptyRows(1000);
        end
        
        function p = loadLabelMat(p, inFileName, scaleFactor)
            tmpLabelMat = imread(inFileName);
            if scaleFactor == 1
                p.dapiMask = tmpLabelMat;
            else
                p.dapiMask = imresize(tmpLabelMat, scaleFactor, 'nearest');
            end
            
            rp = regionprops(p.dapiMask);
            area = [rp.Area];
            idx = area >= p.minNucleusSize;
            rp = rp(idx);
            centroids = [rp.Centroid];
            centroids = round(reshape(centroids,2,[])');
            centroids = single(centroids);
            area = single([rp.Area]);
            
            numCentroids = size(centroids,1);
            status = true(numCentroids, 1);
            maskID = single(zeros(numCentroids, 1));
            colors = single(zeros(numCentroids,3));

            p.nuclei = table((1:numCentroids)',centroids(:,2),centroids(:,1), status, maskID, area', colors,...
                'VariableNames', {'nucID', 'x', 'y', 'status', 'maskID', 'nucleusArea', 'colors'});
            p.addEmptyRows(1000);
        end
        
        function p = stitchCellPoseMasks(p, cpTilePositionsFile, inDir, scaleFactor)
            %Load tiled cellpose outlines from inDir and merge overlapping masks. 
            %Unfortunately, very slow for large scans. Consider writing new function. 
            cpTilePositions = readtable(cpTilePositionsFile, 'ReadRowNames', true);
            rows = regexp(cpTilePositions.tile, '\d+(?=_)', 'Match');
            rows = unique(cellfun(@str2num, [rows{:}]));
            cols = regexp(cpTilePositions.tile, '(?<=_)\d+', 'Match');
            cols = unique(cellfun(@str2num, [cols{:}]));
            tileChannel= regexp(cpTilePositions.tile{1}, '[a-zA-Z]+(?=\d)', 'Match');
            tileChannel = tileChannel{1};
            polyTable = table('Size', [1,4],...
                'VariableTypes', {'single', 'single', 'single', 'polyshape'},...
                'VariableNames', {'polyID', 'row', 'col', 'poly'});
            warning('off', 'MATLAB:polyshape:repairedBySimplify')
            nRows = max(rows);
            nCols = max(cols);
            for i = 1:nRows
                fprintf('Processing row %d of %d\n', i, nRows)
                for ii = 1:nCols
                    if ~mod(ii, 4)
                        fprintf('Processing column %d of %d\n', ii, nCols)
                    end
                    tileName = sprintf('%s%d_%d', tileChannel, rows(i), cols(ii));
                    tilePolys = d2utils.parseCellposeOutlines(fullfile(inDir,sprintf('%s_cp_outlines.txt', tileName)), 'position', cpTilePositions{tileName,{'colStart', 'rowStart'}}, 'scaleFactor', scaleFactor);
                    if ~isempty(tilePolys)
                        polyNeighborIdx = polyTable.row >= i-1 & polyTable.row <= i+1 ... %We may be able to shrink this neighborhood
                            & polyTable.col >= ii-1 & polyTable.col <= ii+1;
                        polyNeighborTable = polyTable(polyNeighborIdx, :);
                        polyOverlapMat = overlaps(polyNeighborTable.poly, tilePolys);
                        if any(polyOverlapMat, 'all')
                            overlapNeighbors = any(polyOverlapMat, 2);
                            tilePolys = regions(union([polyNeighborTable.poly(overlapNeighbors); tilePolys']));
                            tilePolys = tilePolys';
                            polyNeighborIdx(polyNeighborIdx) = overlapNeighbors;
                            polyTable.poly(polyNeighborIdx) = polyshape;
                        end
                        startID = max(polyTable.polyID)+1;
                        tileTable = table((startID:startID+numel(tilePolys)-1)', repmat(i, numel(tilePolys), 1), repmat(ii, numel(tilePolys), 1), tilePolys',...
                            'VariableNames', {'polyID', 'row', 'col', 'poly'});
                        polyTable = [polyTable; tileTable];
                    end
                end
            end
            polyArray = scale(polyTable.poly, scaleFactor);
            [polyX, polyY] = centroid(polyArray);
            polyArea = single(area(polyArray));
            %Remove empty polys
            polyX(polyArea==0) = [];
            polyY(polyArea==0) = [];
            polyArea(polyArea==0) = [];
            
            status = true(numel(polyX), 1);
            maskID = single(zeros(numel(polyX), 1));
            colors = single(zeros(numel(polyX),3));

            p.nuclei = table((1:numel(polyX))',single(polyY),single(polyX), status, maskID, polyArea, colors,...
                'VariableNames', {'nucID', 'x', 'y', 'status', 'maskID', 'nucleusArea', 'colors'});
            p.addEmptyRows(1000);
        end
        
        function p = stitchCellPoseMasks2(p, cpTilePositionsFile, inDir, scaleFactor)
            %Plan of attack. Create empty stitch. Load label matrices in
            %parallel. Cut the perimiter from each mask. Insert tile into stitch matrix. 
            %Find connected components (centroids and area). Add boundary to area if desired. 
            cpTilePositions = readtable(cpTilePositionsFile, 'ReadRowNames', true);
            tiles = cpTilePositions.tile;
            maskArray = cell(1, numel(tiles));
            parfor i = 1:numel(tiles) %Need to test whether this parallelization is faster than just 1 for loop. 
                tileLabelMat = imread(fullfile(inDir, sprintf('%s_cp_masks.tif', tiles{i})));
                tileMask = zeros(size(tileLabelMat));
                cellIDs = unique(tileLabelMat);
                cellIDs(cellIDs == 0) = [];
                for ii = 1:numel(cellIDs)
                    cellMask = tileLabelMat == cellIDs(ii);
                    cellPerim = bwperim(cellMask);
                    cellMask(cellPerim) = 0;
                    tileMask = or(tileMask, cellMask);
                end
                maskArray{i} = imresize(tileMask, scaleFactor);
            end
            maskMatrix = zeros(p.scanObj.stitchDim); %Assuming stitch is the same size as the scan. Alternatively, could use max(rowStart)+height, max(colStart)+width, 
            tileSize = cellfun(@(x) size(x) - [1,1], maskArray, 'UniformOutput', false);
            for i = 1:numel(maskArray)
                left = cpTilePositions.rowStart(i);
                top = cpTilePositions.colStart(i);
                maskMatrix(left:left+tileSize{i}(1),top:top+tileSize{i}(2)) = maskArray{i};
            end
            
            CC = bwconncomp(maskMatrix);
            rp = regionprops(CC);
            area = [rp.Area];
            idx = area >= p.minNucleusSize;
            rp = rp(idx);
            centroids = [rp.Centroid];
            centroids = round(reshape(centroids,2,[])');
            centroids = single(centroids);
            area = single([rp.Area]);
            
            numCentroids = size(centroids,1);
            status = true(numCentroids, 1);
            maskID = single(zeros(numCentroids, 1));
            colors = single(zeros(numCentroids,3));

            p.nuclei = table((1:numCentroids)',centroids(:,2),centroids(:,1), status, maskID, area', colors,...
                'VariableNames', {'nucID', 'x', 'y', 'status', 'maskID', 'nucleusArea', 'colors'});
            p.addEmptyRows(1000);
        end
        
        function p = addEmptyRows(p, n)
            newRows = table('Size', [n, size(p.nuclei, 2)],...
                'VariableNames', p.nuclei.Properties.VariableNames,...
                'VariableTypes', varfun(@class,p.nuclei,'output','cell'));
            newRows.status = false(n, 1);
            newRows.colors = single(zeros(n,3));
            p.nuclei = [p.nuclei; newRows];
        end
        
        function [outNuclei, idx] = getNucleiInRect(p,rect) %rect specified as [x y nrows ncols]
            
            idx = p.nuclei.x >= rect(1) & p.nuclei.x < rect(1) + rect(3) ...
                & p.nuclei.y >= rect(2) & p.nuclei.y < rect(2) + rect(4);
            
            outNuclei = p.nuclei(idx,:);
        end
        
        
        function outNuclei = getValidNucleiInRect(p, rect)
            
            idx = p.nuclei.status ... 
                & p.nuclei.x >= rect(1) & p.nuclei.x < rect(1) + rect(3) ...
                & p.nuclei.y >= rect(2) & p.nuclei.y < rect(2) + rect(4);
            
            outNuclei = p.nuclei(idx,:);
            
        end
        
        function idx = getNucleiInRectIdx(p,rect) %rect specified as [x y nrows ncols]
            
            idx = p.nuclei.x >= rect(1) & p.nuclei.x < rect(1) + rect(3) ...
                & p.nuclei.y >= rect(2) & p.nuclei.y < rect(2) + rect(4) ...
                & p.nuclei.status;
        end
        
        function outNuclei = getNucleiNearRect(p,rect, radius) %rect specified as [x y height width]
            
            idx = p.nuclei.x >= rect(1) - radius & p.nuclei.x < rect(1) + rect(3) + radius ...
                & p.nuclei.y >= rect(2) - radius & p.nuclei.y < rect(2) + rect(4) + radius ...
                & p.nuclei.status;
            
            outNuclei = p.nuclei(idx,:);
        end
        
        function outRect = getRectAroundNuc(p, ID, radius)
            startPos = max([1, 1], p.nuclei{p.nuclei.nucID == ID, {'x', 'y'}} - radius); %Avoid rect outside range of scan. 
            startPos = min(size(p.scanObj.dapiStitch)- radius+1,startPos);
            outRect = [startPos, [2, 2] * radius];
        end
        
        function p = addColors(p)
            randomColors  = single(d2utils.distinguishable_colors(50));
            %randomColors = randomColors(randperm(50), :);
            p.nuclei.colors(p.nuclei.status, :) = [repmat(randomColors, floor(sum(p.nuclei.status)/size(randomColors, 1)), 1);...
                randomColors(1:mod(sum(p.nuclei.status), size(randomColors, 1)), :)];
        end
        
        function p = addCell(p, x, y)
            
            if ~any(p.nuclei.nucID == 0)
                p.addEmptyRows(1000);
            end
            
            tempMaskID = max(p.nuclei.nucID)+1;

            newNuc = table(single(tempMaskID:tempMaskID+numel(x)-1)', single(x), single(y), true(numel(x), 1), single(zeros(numel(x), 1)), single(zeros(numel(x), 1)), single(rand(numel(x), 3)),...
                'VariableNames', p.nuclei.Properties.VariableNames); %Should area be NaN?
            startIdx = find(p.nuclei.nucID == 0, 1, 'first');
            p.nuclei(startIdx:startIdx+numel(x)-1,:) = newNuc;
        end
        
        function p = removeCell(p, x, y)
            if ~isempty(p.nuclei)
                [idx, dist] = knnsearch(p.nuclei{:,{'x','y'}}, [x, y], 'K', 1, 'Distance', 'euclidean');
                if dist < 40 %make sure that the query point is near a nucleus
                    p.nuclei(idx,:) = [];
                end
                
            else
                disp("NucleiTable is empty.")
            end
        end
        
        function p = updateAllMasks(p) %This will overwrite previous maskIDs in nuclei. 
            
            % maskTable = p.maskObj.masks(p.maskObj.masks.dapi,:); % can do
            % this when we rename dapi channel
%             dapiChannelName=p.scanObj.channels{ismember(p.scanObj.channelTypes,'dapi')};
            maskTable = p.maskObj.masks(p.maskObj.masks.dapi,:);
            maskIDs = unique(maskTable.maskID);
            maskIDs(maskIDs == 0) = [];
            p.nuclei = p.nuclei(~(p.nuclei.nucID == 0),:);
            p.nuclei.maskID(:) = single(0);
            p.nuclei.status(:) = true;
            for i = 1:numel(maskIDs)
                idx = inpolygon(p.nuclei.x, p.nuclei.y,...
                    maskTable{maskTable.maskID == maskIDs(i), 'x'}, maskTable{maskTable.maskID == maskIDs(i), 'y'}) & p.nuclei.status;
                p.nuclei.maskID(idx) = maskIDs(i);
                p.nuclei.status(idx) = false;
            end
            p.addEmptyRows(1000);
        end
        
        function p = addNewMask(p, maxCellMask)
%             maxCellMask = max(p.maskObj.masksBB{p.maskObj.masksBB.dapi,'maskID'});
            maskBB = p.maskObj.masksBB{p.maskObj.masksBB.maskID == maxCellMask,{'BB'}}; %Only query nuclei within mask bouding box
            %polyRect = d2utils.boundingCorners2Rect(maskBB);
            
            nucIdx = p.getNucleiInRectIdx(maskBB);
            polyIdx = inpolygon(p.nuclei.x(nucIdx), p.nuclei.y(nucIdx),...
                p.maskObj.masks{p.maskObj.masks.maskID == maxCellMask, 'x'}, p.maskObj.masks{p.maskObj.masks.maskID == maxCellMask, 'y'});
            
            nucIdx(nucIdx) = polyIdx; %Only nuclei in polygon remain true
            p.nuclei.maskID(nucIdx) = maxCellMask;
            p.nuclei.status(nucIdx) = false;
            
        end
       
        function p = updateMasksInRect(p, localRect) 
            %Use for removing masks or if we want to change multiple masks before updating nuclei
            %Probably could use some speeding up. 
            masksInRect = p.maskObj.getChannelMasksInRect(localRect, 'dapi');
            maskIDsinRect = unique(masksInRect.maskID);
            [tmpNuclei, nucIdx] = p.getNucleiInRect(localRect); %It might be faster to not subset the nuclei and just run inpolygon on entire nuclei table. 
            
            %Resest status for nuclei
            tmpNuclei.maskID(:) = single(0);
            tmpNuclei.status(:) = true;
            tmpNuclei.status(tmpNuclei.nucID == 0) = false;
            for i = 1:numel(maskIDsinRect)
                idx = inpolygon(tmpNuclei.x, tmpNuclei.y,...
                    masksInRect{masksInRect.maskID == maskIDsinRect(i), 'x'}, masksInRect{masksInRect.maskID == maskIDsinRect(i), 'y'}) & tmpNuclei.status;
                tmpNuclei.maskID(idx) = maskIDsinRect(i);
                tmpNuclei.status(idx) = false;
                
            end
            p.nuclei.maskID(nucIdx) = tmpNuclei.maskID;
            p.nuclei.status(nucIdx) = tmpNuclei.status;
 
        end
        
        function p = removeMasks(p) 
            
            masksToRemove = setdiff(p.nuclei.maskID, p.maskObj.masksBB.maskID(p.maskObj.masksBB.dapi));
            masksToRemove(masksToRemove == 0) = [];
            if ~isempty(masksToRemove)
                nucIdx = ismember(p.nuclei.maskID, masksToRemove);
                p.nuclei.maskID(nucIdx) = single(0);
                p.nuclei.status(nucIdx) = true;
                p.nucleiChanged = true;
            end
            
        end
        
        function p = removeMasksByID(p, masksToRemove)
            
            nucIdx = ismember(p.nuclei.maskID, masksToRemove);
            if any(nucIdx)
                p.nuclei.maskID(nucIdx) = single(0);
                p.nuclei.status(nucIdx) = true;
                p.nucleiChanged = true;
            end
            
        end
        
        function p = removeMasks2(p, rect)
            %Not sure if it'll be faster to first get nuclei and masks in
            %rect. 
            [nucInRect, nucIdx] = p.getNucleiInRect(rect);
            maskIDsInRect = p.maskObj.getChannelMaskIDsInRect(rect, 'dapi');
            maskIDsInRect(maskIDsInRect == 0) = [];
            goodNucIdx = ~ismember(nucInRect.maskID, maskIDsInRect);
            nucIdx(nucIdx) = goodNucIdx;
            p.nuclei.maskID(nucIdx) = single(0);
            p.nuclei.status(nucIdx) = true;
            p.nucleiChanged = true;
        end
        
        function saveNucleiTable(p)
            if ~isempty(p.nuclei)
                writetable(p.nuclei(~(p.nuclei.nucID == 0),{'nucID', 'x', 'y', 'status', 'maskID', 'nucleusArea'}), p.nucleiFile) %Not saving nuclei.colors
            else
                fprintf("nuclei is empty. Run findNuclei and try again")
            end
        end
        
    end
    
end 