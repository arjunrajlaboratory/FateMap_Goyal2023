classdef cellTable < handle
    
    properties (Access = public)
        
        cells
        maskObj %Mask object 
        minNucleusSize = 1000; % Should set method call findCells? 
        dapiMask

    end
    
    methods
        
        function p = cellTable(maskObj, varargin)
            p.maskObj = maskObj;
            if nargin == 1
                fprintf('New Table\n');
                p.cells = cell2table(cell(0,6), 'VariableNames', {'cellID', 'x', 'y', 'status', 'maskID', 'nucleusArea'}); 
            elseif nargin == 2
                fprintf('Loading Table\n');
                p.cells = readtable(varargin{1},'TextType','string');
            end
        end
        
        function p = stitchDAPImask(p, scanObject, varargin)  
            tileTable = scanObject.tilesTable;
            tilesTmp = transpose(scanObject.scanMatrix);
            tiles = tilesTmp(:);
            [height, width] = scanObject.tileSize();
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'logical');
            channel = find(scanObject.channels == 'dapi');
            reader = bfGetReader(scanObject.scanFile);
            iPlane = reader.getIndex(0, channel - 1, 0) + 1;
            
            if nargin == 2
                s = 0.1;
            elseif nargin == 3
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
            
            p.dapiMask = tmpStitch;
        end
        
        function p = findCells(p)
            
            CC = bwconncomp(p.dapiMask);
            rp = regionprops(CC);
            area = [rp.Area];
            idx = area >= p.minNucleusSize;
            rp = rp(idx);
            centroids = [rp.Centroid];
            centroids = round(reshape(centroids,2,[])');
            centroids = single(centroids);
            
            area = single([rp.Area]);
            
            status = true(height(centroids), 1);
            maskID = single(zeros(height(centroids), 1));

            p.cells = table((1:height(centroids))',centroids(:,2),centroids(:,1), status, maskID, area',...
                'VariableNames', {'cellID', 'x', 'y', 'status', 'maskID', 'nucleusArea'});
        end
        
        function [outCells, idx] = getCellsInRect(p,rect) %rect specified as [x y nrows ncols]
            
            idx = p.cells.x >= rect(1) & p.cells.x < rect(1) + rect(3) ...
                & p.cells.y >= rect(2) & p.cells.y < rect(2) + rect(4);
            
            outCells = p.cells(idx,:);
        end
        
        
        function outCells = getValidCellsInRect(p, rect)
            
            idx = p.cells.status ... 
                & p.cells.x >= rect(1) & p.cells.x < rect(1) + rect(3) ...
                & p.cells.y >= rect(2) & p.cells.y < rect(2) + rect(4);
            
            outCells = p.cells(idx,:);
            
        end
        
        function idx = getCellsInRectIdx(p,rect) %rect specified as [x y nrows ncols]
            
            idx = p.cells.x >= rect(1) & p.cells.x < rect(1) + rect(3) ...
                & p.cells.y >= rect(2) & p.cells.y < rect(2) + rect(4) ...
                & p.cells.status;
        end
        
        function outCells = getCellsNearRect(p,rect, radius) %rect specified as [x y nrows ncols]
            
            idx = p.cells.x >= rect(1) - radius & p.cells.x < rect(1) + rect(3) + radius ...
                & p.cells.y >= rect(2) - radius & p.cells.y < rect(2) + rect(4) + radius ...
                & p.cells.status;
            
            outCells = p.cells(idx,:);
        end
        
        function p = addCell(p, x, y)  
            if ~isempty(p.cells)
                maxID = max(p.cellID);
                newCell = {maxID+1, single(x), single(y), true, single(0), single(0)}; %Should area be NaN?
                p.cells = [p.cells; newCell];
                
            else
                p.cells(1,:) = {single(1), single(x), single(y), true, single(0), single(area)};
            end
        end
        
        function p = removeCell(p, x, y)
            if ~isempty(p.cells)
                [idx, dist] = knnsearch(p.cells{:,{'x','y'}}, [x, y], 'K', 1, 'Distance', 'euclidean');
                if dist < 40 %make sure that the query point is near a cell
                    p.cells(idx,:) = [];
                end
                
            else
                disp("cells is empty.")
            end
        end
        
        function p = updateAllMasks(p) %This will overwrite previous maskIDs in cellTable. 
            
            maskTable = p.maskObj.masks(p.maskObj.masks.dapi);
            maskIDs = unique(maskTable.maskID);
           
            for i = 1:numel(maskIDs)
                idx = inpolygon(p.cells.x, p.cells.y,...
                    maskTable{maskTable.maskID == maskIDs(i), 'x'}, maskTable{maskTable.maskID == maskIDs(i), 'y'});
                p.cells.maskID(idx) = maskIDs(i);
                p.cells.status(idx) = false;
                
            end
            
        end
        
        function p = addNewMask(p)
            
            maxCellMask = max(p.maskObj.masks{p.maskObj.masks.dapi,'maskID'});
            idx = inpolygon(p.cells.x, p.cells.y,...
                    p.maskObj.masks{p.maskObj.masks.maskID == maxCellMask, 'x'}, p.maskObj.masks{p.maskObj.masks.maskID == maxCellMask, 'y'}) & p.cells.status;
            p.cells.maskID(idx) = maxCellMask;
            p.cells.status(idx) = false;
            
        end
       
        function p = updateMasksInRect(p, localRect) %Use for removing masks or if we want to change multiple masks before updating cellTable
            masksInRect = p.maskObj.getChannelMasksInRect(localRect, 'dapi');
            maskIDsinRect = unique(masksInRect.maskID);
            [tmpCells, cellIdx] = p.getCellsInRect(localRect); %It might be faster to not subset the cells and just run inpolygon on entire cellTable. 
            
            %Resest status for cellTable
            p.cells.maskID(cellIdx) = single(0);
            p.cells.status(cellIdx) = true;
            for i = 1:numel(maskIDsinRect)
                idx = inpolygon(tmpCells.x, tmpCells.y,...
                    masksInRect{masksInRect.maskID == maskIDsinRect(i), 'x'}, masksInRect{masksInRect.maskID == maskIDsinRect(i), 'y'}) & tmpCells.status;
                tmpCells.maskID(idx) = maskIDs(i);
                tmpCells.status(idx) = false;
                
            end
            p.cells.maskID(cellIdx) = tmpCells.maskID;
            p.cells.status(cellIdx) = tmpCells.status;
 
        end
        
        function p = removeMasks(p) 
            %If cell falls within multiple masks and only 1 is removed,
            %this function may incorrectly set status to true. Prefer using
            %updateMasksInRect to remove masks
            
            masksToRemove = setdiff(p.cells.maskID, p.maskObj.masks.maskID(p.maskObj.masks.dapi));
            p.cells.maskID(ismember(p.cells.maskID, masksToRemove)) = single(0);
            p.cells.status(ismember(p.cells.maskID, masksToRemove)) = true;
        end
        
        
        function [] = saveCellTable(p, varargin)
            if ~isempty(p.cells)
                if nargin == 1
                    writetable(p.cells, 'cells.csv')
                elseif nargin ==2 
                    writetable(p.cells, varargin{1})
                end
            else
                fprintf("cells is empty. Run findCells and try again")
            end
        end
        
    end
    
end 