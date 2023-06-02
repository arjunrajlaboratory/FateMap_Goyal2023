classdef maskTable < handle
    
    properties (Access = public)
        
        masks
        masksBB
        channels
        heightMaskTable = 100000 %Used to preallocate rows for maskTable
        heightMaskTableBB = 100000
        masksFile = 'masks.csv'
    end
    
    methods
        
        function p = maskTable(scanObject, varargin)
            %channels = convertStringsToChars(scanObject.channels);
%             p.channels = channels; % all channels can have mask
            p.channels = [scanObject.channels(ismember(upper(scanObject.channelTypes),'FISH')), {'dapi'}]; %There can only be masks for spots (FISH channels) and cells (dapi) 
            if nargin == 1
                fprintf('New Mask Table\n');
                p.masks = table('Size', [p.heightMaskTable, numel(p.channels) + 3],...
                    'VariableNames', [{'maskID','x','y'} , p.channels],...
                    'VariableTypes', [repmat({'single'}, 1, 3) , repmat({'logical'}, 1, numel(p.channels))]);
                p.masksBB = table('Size', [p.heightMaskTableBB, numel(p.channels) + 2],...
                    'VariableNames', [{'maskID','BB'} , p.channels],...
                    'VariableTypes', [repmat({'single'}, 1, 2) , repmat({'logical'}, 1, numel(p.channels))]);
                p.masksBB.BB = single(zeros(p.heightMaskTableBB, 4));
%                 p.masks = table(false(0,numel(p.channels) + 3), 'VariableNames', [{'maskID','x','y'} , p.channels]);
%                 p.masksBB = array2table(false(0,numel(p.channels) + 5), 'VariableNames', [{'maskID','x','y', 'h', 'w'} , p.channels]);
            elseif nargin == 2 % Otherwise, load the specified table
                fprintf('Loading Table\n');
                p.masksFile = varargin{1};
                tmpMasks = readtable(varargin{1},'TextType','string');
                tmpMasks = convertvars(tmpMasks,{'maskID', 'x', 'y'},'single');
                p.masks = convertvars(tmpMasks,4:size(tmpMasks, 2),'logical');
                p.allMasks2BB();
                p.addEmptyMasks();
                p.addEmptyMasksBB();
            end
        end
        
        function p = addMask(p,maskPoly,localRect,channel)
            
            tempMaskID = single(max(p.masks.maskID)+1);
            
             if sum(p.masks.maskID == 0) < size(maskPoly, 1)
                p.addEmptyMasks();               
            end
            
            if ~any(p.masksBB.maskID == 0)
                p.addEmptyMasksBB();                
            end
            
            [x,y] = d2utils.localToGlobalCoords(localRect,maskPoly(:,2),maskPoly(:,1));
            
            channelIdx = ismember(p.channels, channel);
            
            tmpCoords = table(repmat(tempMaskID,  length(x), 1), single(x), single(y), 'VariableNames',{'maskID', 'x', 'y'});
            tmpChannelTable = array2table(repmat(channelIdx, length(x), 1), 'VariableNames', p.channels);
            startIdx = find(p.masks.maskID == 0, 1, 'first');

            p.masks(startIdx:startIdx+length(x)-1,:) = [tmpCoords, tmpChannelTable];

            BB = d2utils.polygonBoundingBox([x,y]);
            tmpCoords = cell2table({tempMaskID, single(BB)}, 'VariableNames',{'maskID', 'BB'});
            startIdx = find(p.masksBB.maskID == 0, 1, 'first');
            p.masksBB(startIdx,:) = [tmpCoords, tmpChannelTable(1,:)];

        end
        
        function tempMaskID = addMaskLocalCoords(p,maskPoly,channelOrChannels)
            
            tempMaskID = single(max(p.masks.maskID)+1);
            
            numMasks = size(maskPoly,1);
            if sum(p.masks.maskID == 0) < numMasks
                p.addEmptyMasks();                
            end
            
            if ~any(p.masksBB.maskID == 0)
                p.addEmptyMasksBB();                
            end
            
%             [x,y] = d2utils.localToGlobalCoords(localRect,maskPoly(:,2),maskPoly(:,1));
            
            channelIdx = ismember(p.channels, channelOrChannels);
            
            tmpCoords = table(repmat(tempMaskID,  numMasks, 1), single(maskPoly(:,2)), single(maskPoly(:,1)), 'VariableNames',{'maskID', 'x', 'y'});
            tmpChannelTable = array2table(repmat(channelIdx, numMasks, 1), 'VariableNames', p.channels);
            startIdx = find(p.masks.maskID == 0, 1, 'first');
            p.masks(startIdx:startIdx+numMasks-1,:) = [tmpCoords, tmpChannelTable];
            
            BB = d2utils.polygonBoundingBox(fliplr(maskPoly));
            tmpCoords = cell2table({tempMaskID, single(BB)}, 'VariableNames',{'maskID', 'BB'});
            startIdx = find(p.masksBB.maskID == 0, 1, 'first');
            p.masksBB(startIdx,:) = [tmpCoords, tmpChannelTable(1,:)];

        end
        
        function p = addEmptyMasks(p)
            newRows = table('Size', [p.heightMaskTable, numel(p.channels) + 3],...
                'VariableNames', p.masks.Properties.VariableNames,...
                'VariableTypes', varfun(@class,p.masks,'output','cell'));
            p.masks = [p.masks; newRows]; 
        end
        
        function p = addEmptyMasksBB(p)
            newRows = table('Size', [p.heightMaskTableBB, numel(p.channels) + 2],...
                'VariableNames', p.masksBB.Properties.VariableNames,...
                'VariableTypes',  varfun(@class,p.masksBB,'output','cell'));
            newRows.BB = zeros(p.heightMaskTableBB, 4);
            p.masksBB = [p.masksBB; newRows]; 
        end
        
        function p = removeMasks(p,maskIDs)
             p.masks(ismember(p.masks.maskID,maskIDs),:) = [];
             p.masksBB(ismember(p.masksBB.maskID,maskIDs),:) = [];
             
        end
        
        function p = updateMaskPoly(p,channel,maskID,maskPoly,localRect)
            p.removeMasks(maskID);
            p.addMask(maskPoly,localRect,channel)
        end
        
        function outMasks = getAllMasksInRect(p, rect)
            
            idx = rectint(p.masksBB.BB,rect) > 0;
                        
            outMasks = p.masks(ismember(p.masks.maskID, p.masksBB.maskID(idx)) ,:);
        end
        
        function outMasks = getChannelMasksInRect(p, rect, channel)
            
            idx = rectint(p.masksBB.BB, rect) > 0 & p.masksBB{:, channel};
           
            outMasks = p.masks(ismember(p.masks.maskID, p.masksBB.maskID(idx)) ,:);
        end
        
        function outMasks = getDapiAndChannelMasksInRect(p, rect, channel)
            
            idx = rectint(p.masksBB.BB, rect) > 0 & (p.masksBB{:, channel} | p.masksBB{:, 'dapi'});
           
            outMasks = p.masks(ismember(p.masks.maskID, p.masksBB.maskID(idx)) ,:);
        end
                
        function outMaskIDs = getChannelMaskIDsInRect(p, rect, channel)
            
            idx = rectint(p.masksBB.BB, rect) > 0 & p.masksBB{:, channel};
           
            outMaskIDs = p.masksBB.maskID(idx);
        end
        
        function p = removeMasksByPoints(p, points, rect, channel)
            masksInRect = p.getChannelMasksInRect(rect, channel);
            maskIDs = unique(masksInRect.maskID);
            maskIDs(maskIDs==0) = [];
            [x,y] = d2utils.localToGlobalCoords(rect,points(:,2),points(:,1));
            for i = 1:numel(maskIDs)
                if any(inpolygon(x, y, masksInRect{masksInRect.maskID == maskIDs(i), 'x'}, masksInRect{masksInRect.maskID == maskIDs(i), 'y'}))
                    p.removeMasks(maskIDs(i));
                end
            end
        end
        
        function varargout = removeMasksByLocalPoints(p, points, rect)
            masksInRect = p.getAllMasksInRect(rect);
            maskIDs = unique(masksInRect.maskID);
            maskIDs(maskIDs==0) = [];
            removedMaskIDs = [];
            globalMask = [];
            removedMaskBB = {};
            for i = 1:numel(maskIDs)
                if any(inpolygon(points(:,2), points(:,1), masksInRect{masksInRect.maskID == maskIDs(i), 'x'}, masksInRect{masksInRect.maskID == maskIDs(i), 'y'}))
                    removedMaskIDs = [removedMaskIDs, maskIDs(i)];
                    tmpMaskBB = p.masksBB(p.masksBB.maskID == maskIDs(i),:);
                    removedMaskBB = [removedMaskBB; tmpMaskBB{:,'BB'}];
                    globalMask = [globalMask, sum(tmpMaskBB{:,3:end}) > 1]; %There may be a better way of doing this
                    p.removeMasks(maskIDs(i));
                end
            end
            varargout{1} = removedMaskIDs;
            varargout{2} = globalMask;
            varargout{3} = removedMaskBB;
        end
        
        function p = allMasks2BB(p)
            p.masksBB = table('Size', [0, numel(p.channels) + 2],...
                    'VariableNames', [{'maskID','BB'} , p.channels],...
                    'VariableTypes', [repmat({'single'}, 1, 2) , repmat({'logical'}, 1, numel(p.channels))]);
            uniqueMaskIDs = unique(p.masks.maskID);
            uniqueMaskIDs(uniqueMaskIDs==0) = [];
            for i = 1:numel(uniqueMaskIDs)
                tempMask = p.masks(p.masks.maskID == uniqueMaskIDs(i), :);
                BB = d2utils.polygonBoundingBox(tempMask{:,{'x', 'y'}});
                tmpCoords = cell2table({uniqueMaskIDs(i), single(BB)}, 'VariableNames',{'maskID', 'BB'});
                p.masksBB = [p.masksBB; [tmpCoords, tempMask(1,4:end)]];
            end
        end
       
        function saveMasksTable(p)
           if ~isempty(p.masks)
               writetable(p.masks(~(p.masks.maskID) == 0,:), p.masksFile);
           else
               fprintf("Masks table is empty. Not saving masks.csv")
           end
        end
        
    end
    
end 