classdef spotTable < handle
    
    properties (Access = public)
        
        scanObj
        maskObj 
        nucleiObj 
        
        spots
        centroidLists %Not sure if this should go into the View
        
        sigma = 2;
        thresholds = [];
        minChannelIntensities
        minChannelIntensityIsBlockSpecific
        spotChannels
        maxDistance = 200;  
        theFilter
        percentileToKeep = 98;
        threshFactor = 2;
        spotsIntensitiesWithMasked
        spotsIntensitiesNoMasked
        expressionColorPal = {'BuYlRd', 'YlOrRd', 'GrBu', 'BuGnYlRd'}
        paletteIdx = 1;
        borderSpotIdx;
        spotsFile %Output filenames
        spotsSummaryFile = 'spotsSummary.csv'
    end
    
    methods
        
        function p = spotTable(scanObject, maskObject, nucleiObject, varargin)

            n=inputParser;
            n.addRequired('scanObject', @(x) mustBeA(x, 'scanObject'))
            n.addRequired('maskObject', @(x) mustBeA(x, 'maskTable'))
            n.addRequired('nucleiObject', @(x) mustBeA(x, 'nucleiTable'))
            n.addParameter('spotsFile','',@ischar)
            n.addParameter('sigma',[],@(x) validateattributes(x,{'numeric'},{'positive','nonzero'}))
            n.addParameter('thresholds',[], @(x) isempty(x) || (isnumeric(x) && size(x,1)==1) )
            n.addParameter('aTrousMinThreshFactor', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0))
            n.parse(scanObject, maskObject, nucleiObject, varargin{:});
            
            p.scanObj = scanObject;
            p.maskObj = maskObject;
            p.nucleiObj = nucleiObject;
            
            p.parseScanSummary(p.scanObj.scanSummaryFile);

            if ~isempty(n.Results.thresholds)
                p.userInputThresholds(n.Results.thresholds) % will set p.thresholds. Otherwise p.thresholds might still be nonempty if it came from scanSummary.txt
            end

            if ~isfile(n.Results.spotsFile)
                fprintf('Unable to find %s in your current directory. Creating a new spots object\n', n.Results.spotsFile)
%                 fprintf('New spots table\n');
                p.spotsFile = n.Results.spotsFile;
                % even if we have params sigma, thresholds,
                % aTrousMinThreshFactor in scanSummary.txt, overwrite with user inputs
                if ~isempty(n.Results.sigma)
                    p.sigma=n.Results.sigma;
                end
                
                if ~isempty(n.Results.aTrousMinThreshFactor)
                    p.threshFactor=n.Results.aTrousMinThreshFactor;
                end

                if isempty(p.thresholds)
                    p.minChannelIntensityIsBlockSpecific=true; % findSpots5 will return spots above minIntensity=autoThresholdForThisBlock/p.threshFactor
                    p.minChannelIntensities=[]; % don't really need to set this, already the case by default
                else
                    p.minChannelIntensityIsBlockSpecific=false; % findSpots5 will return spots above minIntensity=globalThreshold/p.threshFactor (where globalThreshold is the channel-specific value from p.thresholds, which was loaded from scanSummary.txt or was user-input
                    p.minChannelIntensities=p.thresholds/p.threshFactor;
                end
                %Commenting out the code below in order to more easily create empty spotTable objects.  
%                 disp('Finding spots. This may take several minutes.')
%                 p.findSpots5(); %Run before contrasting scans
%                 p.maskBorderSpots();
%                 disp('Finished finding spots')
%                 p.assignSpotsToNuclei();
                
            else % there is a spots file here (usually spots.csv)
                fprintf('Loading spot table\n');
                p.spotsFile = n.Results.spotsFile;
                opts = detectImportOptions(n.Results.spotsFile);
                opts = setvartype(opts, {'spotID', 'x', 'y', 'nearestNucID', 'maskID', 'distanceToNuc'}, 'single');
                opts = setvartype(opts, 'channel', 'string');
                opts = setvartype(opts, {'intensity', 'status'}, 'uint16'); %For some reason, when I set 'status' to 'logical' they all go to false. So doing this instead
                p.spots = readtable(p.spotsFile, opts);
                p.spots.status = logical(p.spots.status);
                [spotIdx, nucIdx] = ismember(p.spots.nearestNucID, p.nucleiObj.nuclei.nucID);
                nucIdx(nucIdx == 0) = [];
                p.spots{spotIdx, 'colors'} = p.nucleiObj.nuclei.colors(nucIdx,:);
                
            end
        end
        
        function p = assignSpotsToNuclei(p)
            %At the moment, assigning spots to nearest valid nucleus (i.e. not masked). 
            %If we want to assign spots to nearest ncleus even if that
            %nucleus is masked (status = false). Then change below to 
            %validNuclei = p.nucleiObj.nuclei(~p.nucleiObj.nuclei.nucID == 0,:)
            %below to validNuclei = p.nucleiObj.nuclei(p.nucleiObj.nuclei.status, :);

            validNuclei = p.nucleiObj.nuclei(p.nucleiObj.nuclei.status, :); %There may be a more efficient way fo doing this using only index arrays and not creating a new table. 
            
            [idx, dist] = knnsearch(validNuclei{:,{'x', 'y'}}, p.spots{:,{'x', 'y'}}, 'K', 1, 'Distance', 'euclidean');
            
            p.spots.nearestNucID = validNuclei.nucID(idx);
            p.spots.distanceToNuc = single(dist);
            p.spots.colors = validNuclei.colors(idx, :);
            
        end
        
        function p = assignSpotsToNucleiChannel(p, channel)
            %At the moment, assigning spots to nearest valid nucleus (i.e. not masked). 
            %If we want to assign spots to nearest ncleus even if that
            %nucleus is masked (status = false). Then change below to 
            %validNuclei = p.nucleiObj.nuclei(~p.nucleiObj.nuclei.nucID == 0,:)
            %below to validNuclei = p.nucleiObj.nuclei(p.nucleiObj.nuclei.status, :);

            validNuclei = p.nucleiObj.nuclei(p.nucleiObj.nuclei.status, :); %There may be a more efficient way fo doing this using only index arrays and not creating a new table. 
            spotIdx = p.spots.channel == channel;
            [nucIdx, dist] = knnsearch(validNuclei{:,{'x', 'y'}}, p.spots{spotIdx, {'x', 'y'}}, 'K', 1, 'Distance', 'euclidean');
            
            p.spots{spotIdx,'nearestNucID'} = validNuclei.nucID(nucIdx);
            p.spots{spotIdx, 'distanceToNuc'} = single(dist);
            p.spots{spotIdx , 'colors'} = validNuclei.colors(nucIdx, :);
        end
        
        function p = assignSpotsInRect(p, rect) %Maybe useful for reassigning spots after add/removing cells 
            
            [spotsInRect, spotIdx] = p.getAllSpotsInRect(rect);
            nucleiNearRect = p.nucleiObj.getNucleiNearRect(rect, p.maxDistance);
            if ~isempty(nucleiNearRect) && any(spotIdx)
                [nucIdx, dist] = knnsearch(nucleiNearRect{:,{'x', 'y'}}, spotsInRect{:,{'x','y'}}, 'K', 1, 'Distance', 'euclidean');
                p.spots.nearestNucID(spotIdx) = nucleiNearRect.nucID(nucIdx);
                p.spots.distanceToNuc(spotIdx) = single(dist);
                p.spots.colors(spotIdx, :) = nucleiNearRect.colors(nucIdx, :);
                p.allIntensities(); %Because there could be new valid spots to count.
            end
        end
        
        function idx = getAllSpotsInRectIndex(p,rect) %rect specified as [x y nrows ncols]
            
            idx = p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
        end
        
        function [outSpots,idx] = getAllSpotsInRect(p,rect) %rect specified as [x y nrows ncols]
            
            idx = p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
            outSpots = p.spots(idx,:);
        end
        
        
        function [outSpots, idx] = getValidSpotsInRect(p,channel,rect) %rect specified as [x y nrows ncols]

            idx = p.spots.channel == channel & p.spots.status ... 
                & p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
            outSpots = p.spots(idx,:);
        end
        
        function [outSpots,idx] = getSpotsInRect(p,channel,rect) %rect specified as [x y nrows ncols]
            %outSpots = p.spots(p.spots.channel == channel,:);
            %outSpots = outSpots(outSpots.status,:);

            idx = p.spots.channel == channel ... 
                & p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
            outSpots = p.spots(idx,:);
        end
        
        function idx = getValidSpotsInRectIndex(p,channel,rect) 
            
            idx = p.spots.channel == channel & p.spots.status ... 
                & p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
        end
        
        function idx = getSpotsInRectIndex(p,channel,rect)

            idx = p.spots.channel == channel ... 
                & p.spots.x >= rect(1) & p.spots.x < rect(1) + rect(3) ...
                & p.spots.y >= rect(2) & p.spots.y < rect(2) + rect(4);
            
        end
                
        function p = updateAllMasks(p)  %Note, this will overwrite previous maskIDs in spotTable. 
            for i = 1:numel(p.spotChannels)
                channel = p.spotChannels{i};
                p = updateChannelMasks(p, channel);
            end 
        end
        
        function p = updateChannelMasks(p, channel)  %Note, this will overwrite previous maskIDs in spotTable. 
            
                channelIdx = p.maskObj.masks{:, channel};
                maskTable = p.maskObj.masks(channelIdx,:);
                maskIDs = unique(maskTable.maskID);
                maskIDs(maskIDs == 0) = [];
                spotIdx = p.spots.channel == channel;
                p.spots.maskID(spotIdx) =  single(0);
                p.updateSpotStatus(channel);
                
                validSpotIdx = p.spots.channel == channel & p.spots.status;
                spotTable = p.spots(validSpotIdx,:);
                
                for i = 1:numel(maskIDs)
                    idx = inpolygon(spotTable.x, spotTable.y,...
                        maskTable{maskTable.maskID == maskIDs(i), 'x'}, maskTable{maskTable.maskID == maskIDs(i), 'y'}) & spotTable.status;
                    spotTable.maskID(idx) = maskIDs(i); 
                    spotTable.status(idx) = false; 
                end
                p.spots.maskID(validSpotIdx) = spotTable.maskID;
                p.spots.status(validSpotIdx) = spotTable.status;
        end
        
        function p = addNewMask(p, channel, maxSpotMask)
%             maxSpotMask = max(p.maskObj.masks{p.maskObj.masks{:,channel},'maskID'});
            maskBB = p.maskObj.masksBB{p.maskObj.masksBB.maskID == maxSpotMask,{'BB'}}; %Only query spots within mask bouding box
            %polyRect = d2utils.boundingCorners2Rect(maskBB);
            spotIdx = p.getSpotsInRectIndex(channel, maskBB);
            idx = inpolygon(p.spots.x(spotIdx), p.spots.y(spotIdx),...
                p.maskObj.masks{p.maskObj.masks.maskID == maxSpotMask, 'x'}, p.maskObj.masks{p.maskObj.masks.maskID == maxSpotMask, 'y'});
            
            spotIdx(spotIdx) = idx; %Only spots in polygon remain true
                
            p.spots.maskID(spotIdx) = maxSpotMask;
            p.spots.status(spotIdx) = false;
            %p.spotsIntensitiesNoMasked{ismember(p.spotChannels,channel)} = sort(uint16(p.getIntensitiesNoMasked(channel)));
        end
        
        function p = addNewMaskAllChannels(p, maxSpotMask)
            maskBB = p.maskObj.masksBB{p.maskObj.masksBB.maskID == maxSpotMask,{'BB'}}; %Only query spots within mask bouding box

            spotIdx = p.getAllSpotsInRectIndex(maskBB);
            idx = inpolygon(p.spots.x(spotIdx), p.spots.y(spotIdx),...
                p.maskObj.masks{p.maskObj.masks.maskID == maxSpotMask, 'x'}, p.maskObj.masks{p.maskObj.masks.maskID == maxSpotMask, 'y'});
            
            spotIdx(spotIdx) = idx; %Only spots in polygon remain true
                
            p.spots.maskID(spotIdx) = maxSpotMask;
            p.spots.status(spotIdx) = false;
        end
        
        function p = updateMasksInRect(p, channel, localRect) 
            %Use for removing masks or if we want to change multiple masks
            %before updating spots
            %Probably could use some speeding up. 
            masksInRect = p.maskObj.getChannelMasksInRect(localRect, channel);
            maskIDsinRect = unique(masksInRect.maskID);
            [tmpSpots, spotIdx] = p.getSpotsInRect(channel, localRect);            
            
            %Resest status for nuclei
            tmpSpots.maskID(:) = single(0);
            tmpSpots.status(:) = true;
            for i = 1:numel(maskIDsinRect)
                idxPoly = inpolygon(tmpSpots.x, tmpSpots.y,...
                    masksInRect{masksInRect.maskID == maskIDsinRect(i), 'x'}, masksInRect{masksInRect.maskID == maskIDsinRect(i), 'y'}) & tmpSpots.status;
                
                tmpSpots.maskID(idxPoly) = maskIDsinRect(i);
                tmpSpots.status(idxPoly) = false;
                
            end
            
            p.spots.maskID(spotIdx) = tmpSpots.maskID;
            p.updateSpotStatus(channel);
            %p.spotsIntensitiesNoMasked{ismember(p.spotChannels,channel)} = sort(uint16(p.getIntensitiesNoMasked(channel)));
        end
        
        function p = removeMasks(p, channel) 
            masksToRemove = setdiff(p.spots.maskID, p.maskObj.masks.maskID(p.maskObj.masks{:,channel}));
            p.spots.maskID(ismember(p.spots.maskID, masksToRemove)) = single(0);
            %p.updateSpotStatus(channel);
            %p.spotsIntensitiesNoMasked{ismember(p.spotChannels,channel)} = sort(uint16(p.getIntensitiesNoMasked(channel)));
        end
        
        function p = removeMasksByID(p, masksToRemove)
            p.spots.maskID(ismember(p.spots.maskID, masksToRemove)) = single(0);
            %p.updateSpotStatus(channel);
            %p.spotsIntensitiesNoMasked{ismember(p.spotChannels,channel)} = sort(uint16(p.getIntensitiesNoMasked(channel)));
        end
        
        function p = removeMasks2(p, channel, rect)
            %Not sure if it'll be faster to first get spots and masks in
            %rect. 
            [spotsInRect, spotIdx] = p.getSpotsInRect(channel, rect);
            maskIDsInRect = p.maskObj.getChannelMaskIDsInRect(rect, channel);
            maskIDsInRect(maskIDsInRect == 0) = [];
            goodSpotIdx = ~ismember(spotsInRect.maskID, maskIDsInRect);
            spotIdx(spotIdx) = goodSpotIdx;
            p.spots.maskID(spotIdx) = single(0);
            %p.updateSpotStatus(channel); %For now, done by controller.
            %p.spotsIntensitiesNoMasked{ismember(p.spotChannels,channel)} = sort(uint16(p.getIntensitiesNoMasked(channel)));
        end
        
        function p = removeMasks3(p, rect)
            %Not sure if it'll be faster to first get spots and masks in
            %rect. 
            [spotsInRect, spotIdx] = p.getSpotsInRect(channel, rect);
            maskIDsInRect = p.maskObj.getChannelMaskIDsInRect(rect, channel);
            maskIDsInRect(maskIDsInRect == 0) = [];
            goodSpotIdx = ~ismember(spotsInRect.maskID, maskIDsInRect);
            spotIdx(spotIdx) = goodSpotIdx;
            p.spots.maskID(spotIdx) = single(0);
            %p.updateSpotStatus(channel); %For now, done by controller.
            %p.spotsIntensitiesNoMasked{ismember(p.spotChannels,channel)} = sort(uint16(p.getIntensitiesNoMasked(channel)));
        end
        
        function p = setThreshold(p, channel, value)
            p.thresholds(ismember(p.spotChannels, channel)) = max(1, value); %Threshold >= 1. 
            %Update spot status
            p.updateSpotStatus(channel);
            p.updateCentroidList(channel);
        end
       
        function p = updateSpotStatus(p,channel)
            threshold = p.thresholds(ismember(p.spotChannels, channel));
            channelIdx = ismember(p.spots.channel, channel);
            spotIdx = p.spots.intensity(channelIdx) >= threshold ...
                & p.spots.distanceToNuc(channelIdx) <= p.maxDistance & p.spots.maskID(channelIdx) == 0 ...
                & ismember(p.spots.nearestNucID(channelIdx), p.nucleiObj.nuclei.nucID(p.nucleiObj.nuclei.status));%spots assigned to masked cells will be set to false. 
            
            p.spots.status(channelIdx) = spotIdx;
            %p.spots.status(~spotIdx) = false;
        end
        
        function p = updateAllSpotStatus(p)
            for i = 1:numel(p.spotChannels)
                p.updateSpotStatus(p.spotChannels{i});
            end
        end
       
        function p = updateSpotStatusInRect(p, channel, rect)
            %Add buffer to rect based on max distance allowed from spot to
            %nucleus
            [spotsInRect, spotIdx] = p.getSpotsInRect(channel, rect);
            threshold = p.thresholds(ismember(p.spotChannels, channel));

            goodSpotIdx = spotsInRect.intensity >= threshold ...
                & spotsInRect.distanceToNuc <= p.maxDistance & spotsInRect.maskID == 0 ...
                & ismember(spotsInRect.nearestNucID, p.nucleiObj.nuclei.nucID(p.nucleiObj.nuclei.status));%spots assigned to masked cells will be set to false. 
            
            p.spots.status(spotIdx) = goodSpotIdx;
        end
        
        function p = updateAllSpotStatusInRect(p, rect)
            for i = 1:numel(p.spotChannels)
                p.updateSpotStatusInRect(p.spotChannels{i}, rect);
            end
        end
        
        function p = findSpots3(p) 
            %Find spots on stitched image.
            %Do not run on auto-contrasted stitches. 
            %Only run on non-contrasted stitches. 
            if isempty(p.theFilter) == 2
                filt = -fspecial('log',20,2);
                p.theFilter = filt;
            else
                filt = p.theFilter;
            end
            
            x = [];
            y = [];
            intensity = [];
            channel = [];
            threshPercentile = p.percentileToKeep;
            [splitMats, startCoords] = p.scanObj.splitStitch(p.spotChannels);
            
            
            for i = 1:numel(p.spotChannels)
                fprintf('Finding %s spots\n',p.spotChannels{i});
                splitMat = splitMats(:,i);
                tempX = cell(numel(splitMat), 1);
                tempY = cell(numel(splitMat), 1);
                tempIntensity = cell(numel(splitMat), 1);
                
                
                thresh=p.thresholds(i);
                %parpool('threads')
                parfor ii = 1:numel(splitMat)
                    %Consider changing connectivity to 4 to find more spots in dense areas. 
                    [tempX{ii}, tempY{ii}, tempIntensity{ii}] = d2utils.findSpotsInImage(splitMat{ii}, threshPercentile, 'filter', filt, 'shift', startCoords(ii,:));
                end
                x = [x ; cell2mat(tempX)];
                y = [y ; cell2mat(tempY)];
                intensity = [intensity ; cell2mat(tempIntensity)];
                channel = [channel ; repmat(string(p.spotChannels{i}),size(cell2mat(tempX), 1),1)];
            end 
            
            spotID = single((1:length(x)))';
            nearestNucID = single(zeros(length(x),1));
            maskID = single(zeros(length(x),1));
            status = true(length(x),1);
            dist =  single(zeros(length(x),1));
            p.spots = table(spotID, single(x), single(y), intensity, nearestNucID, status, maskID, channel, dist,...
                'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestNucID', 'status', 'maskID', 'channel', 'distanceToNuc'});
        end
        
        function p = findSpots4(p, varargin) 
            %Find spots on stitched image.
            %Do not run on auto-contrasted stitches. 
            %Only run on non-contrasted stitches. 
            if nargin == 2
                p.threshFactor = varargin{1}; %Global threshFactor set when launching GUI.
            end
            
            x = [];
            y = [];
            intensity = [];
            channel = [];
            
            spotThreshFactor = p.threshFactor; %Assign variables here to reduce overhead in parfor loop below 
            [splitMats, startCoords] = p.scanObj.splitStitch(p.spotChannels);
            for i = 1:numel(p.spotChannels)
                fprintf('Finding %s spots\n',p.spotChannels{i});
                splitMat = splitMats(:,i);
                tempX = cell(numel(splitMat), 1);
                tempY = cell(numel(splitMat), 1);
                tempIntensity = cell(numel(splitMat), 1);
                
                %parpool('threads')
                parfor ii = 1:numel(splitMat)
                    %Consider changing connectivity to 4 to find more spots in dense areas. 
                    [tempX{ii}, tempY{ii}, tempIntensity{ii}] = d2utils.findSpotsaTrous(splitMat{ii}, 'shift', startCoords(ii,:) - 1, 'threshFactor', spotThreshFactor); % the - 1 was added 6-Sep-2021

                end
                x = [x ; cell2mat(tempX)];
                y = [y ; cell2mat(tempY)];
                intensity = [intensity ; cell2mat(tempIntensity)];
                channel = [channel ; repmat(string(p.spotChannels{i}),size(cell2mat(tempX), 1),1)];
            end 
            
            n = length(x);
            spotID = single((1:n))';
            nearestNucID = single(zeros(n,1));
            maskID = single(zeros(n,1));
            status = true(n,1);
            dist =  single(zeros(n,1));
            intensity = uint16(intensity); %Will cut off values <0 and >65535 but very unlikely.  
            p.spots = table(spotID, single(x), single(y), intensity, nearestNucID, status, maskID, channel, dist,...
                'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestNucID', 'status', 'maskID', 'channel', 'distanceToNuc'});
        end
        
        function p = findSpots5(p)
            %Find spots on stitched image.
            %Do not run on auto-contrasted stitches.
            %Only run on non-contrasted stitches.
            
            x = [];
            y = [];
            intensity = [];
            channel = [];
            
            spotThreshFactor = p.threshFactor; %Assign variables here to reduce overhead in parfor loop below 
            spotSigma = p.sigma;
            [splitMats, startCoords] = p.scanObj.splitStitch(p.spotChannels);
            
            if ~isempty(p.thresholds)
                assert(length(p.thresholds)==length(p.spotChannels))
            end
            
            if ~p.minChannelIntensityIsBlockSpecific
                fprintf('spot finder will use these minChannelIntensities for the entire scan:\n\n')
                strlen=max(cellfun(@(x) length(x),p.spotChannels));
                fprintf('%*s   minChannelIntensity\n',strlen,'spotChannel')
                D=[num2cell(repmat(strlen,1,length(p.spotChannels))); reshape(p.spotChannels,1, []); num2cell(reshape(p.minChannelIntensities,1, []))];
                fprintf('%*s     %.2f\n',D{:})
            end
            for i = 1:numel(p.spotChannels)
                fprintf('Finding %s spots\n',p.spotChannels{i});
                splitMat = splitMats(:,i);
                tempX = cell(numel(splitMat), 1);
                tempY = cell(numel(splitMat), 1);
                tempIntensity = cell(numel(splitMat), 1);
                
                if p.minChannelIntensityIsBlockSpecific
                    minRegionalMax=[]; % if minRegionalMax is empty, d2utils.findSpotsaTrous2 then will make minRegionalMax to be block-specific, based on the auto-threshold found for that block
                else
                    minRegionalMax=p.minChannelIntensities(i);
                end
                
                %parpool('threads')
                parfor ii = 1:numel(splitMat)
                    %Consider changing connectivity to 4 to find more spots in dense areas.
                    [tempX{ii}, tempY{ii}, tempIntensity{ii}] = d2utils.findSpotsaTrous2(splitMat{ii},'minRegionalMax',minRegionalMax,'sigma', spotSigma, 'shift', startCoords(ii,:) - 1,'threshFactor', spotThreshFactor); % the - 1 was added 6-Sep-2021
                end
                x = [x ; cell2mat(tempX)];
                y = [y ; cell2mat(tempY)];
                intensity = [intensity ; cell2mat(tempIntensity)];
                channel = [channel ; repmat(string(p.spotChannels{i}),size(cell2mat(tempX), 1),1)];
            end
            
            n = length(x);
            spotID = single((1:n))';
            nearestNucID = single(zeros(n,1));
            maskID = single(zeros(n,1));
            status = true(n,1);
            dist =  single(zeros(n,1));
            intensity = uint16(intensity); %Will cut off values <0 and >65535 but very unlikely.
            p.spots = table(spotID, single(x), single(y), intensity, nearestNucID, status, maskID, channel, dist,...
                'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestNucID', 'status', 'maskID', 'channel', 'distanceToNuc'});
        end
                
        function p = findSpotsChannel(p, channel, varargin) % this function is not updated to call with the recent 
            if nargin == 3
                tmpThreshFactor = varargin{1}; %Channel specific threshFactor
            else
                tmpThreshFactor = p.threshFactor; 
            end
            %Use to refind spots with new threshold.
            %Do not run on auto-contrasted stitches. 
            %Only run on non-contrasted stitches. 
            [splitMats, startCoords] = p.scanObj.splitStitch(channel);
            fprintf('Finding %s spots\n',channel);
            %Below, setting block size to 1000 pixels. Can make this larger if you have very even illumination/spot intensities. Just beware of number of bits.
            splitMat = splitMats(:,1);
            tempX = cell(numel(splitMat), 1);
            tempY = cell(numel(splitMat), 1);
            tempIntensity = cell(numel(splitMat), 1);
            %parpool('threads')
            parfor ii = 1:numel(splitMat)
                %Consider changing connectivity to 4 to find more spots in dense areas.
                [tempX{ii}, tempY{ii}, tempIntensity{ii}] = d2utils.findSpotsaTrous(splitMat{ii}, 'shift', startCoords(ii,:) - 1, 'threshFactor', tmpThreshFactor); % the - 1 was added in Sep 2021
            end
            x = cell2mat(tempX);
            y = cell2mat(tempY);
            intensity = cell2mat(tempIntensity);
            n = size(x, 1);
            channelMat = repmat(channel,n,1);
            
            spotID = single((1:n))';
            nearestNucID = single(zeros(n,1));
            maskID = single(zeros(n,1));
            status = true(n,1);
            dist =  single(zeros(n,1));
            colors = single(zeros(n,3));
            newSpots = table(spotID, single(x), single(y), intensity, nearestNucID, status, maskID, channelMat, dist, colors,...
                'VariableNames', {'spotID', 'x', 'y', 'intensity', 'nearestNucID', 'status', 'maskID', 'channel', 'distanceToNuc', 'colors'});
            %Replace spots for channel
            p.spots(p.spots.channel == channel, :) = [];
            p.spots = [p.spots; newSpots];
        end
        
        function maskBorderSpots(p)
            %Only run this with stitchedScans available
            mask = p.scanObj.stitchedScans.stitches{1} == 0;	
            clearBorder = imclearborder(mask);
            mask = xor(mask, clearBorder);
            maskDilated = imdilate(mask, true(17)); %Adding 8 pixels to borderMask. ~3 pixels are needed for findSpots3. 
            [maskXidx, maskYidx] = find(maskDilated);
            p.borderSpotIdx = ismember(p.spots{:,{'x', 'y'}},[maskXidx, maskYidx], 'rows'); 
            
            p.spots.status(p.borderSpotIdx) = false;
            p.spots.intensity(p.borderSpotIdx) = single(0);%So that the spots aren't unmasked on update. May want an additional flag like maskID = NaN.    
        end
        
        function intensities = getIntensities(p, channel)
            intensities = p.spots{p.spots.channel == channel & p.spots.distanceToNuc <= p.maxDistance, {'intensity'}};
            intensities(intensities == 0) = []; %Removes border spots. 
        end
        
        function intensities = getIntensitiesNoMasked(p, channel)
            intensities = p.spots{p.spots.channel == channel & p.spots.distanceToNuc <= p.maxDistance & p.spots.maskID == 0, {'intensity'}};
            intensities(intensities == 0) = []; %Removes border spots. 
        end
        
        function p = allIntensities(p)
            p.spotsIntensitiesWithMasked = cell(0, numel(p.spotChannels));
            for i = 1:numel(p.spotChannels)
                p.spotsIntensitiesWithMasked{i} = sort(uint16(p.getIntensities(p.spotChannels{i})));
            end
        end
        
        function p = allIntensitiesNoMasked(p)
            p.spotsIntensitiesNoMasked = cell(0, numel(p.spotChannels));
            for i = 1:numel(p.spotChannels)
                p.spotsIntensitiesNoMasked{i} = sort(uint16(p.getIntensitiesNoMasked(p.spotChannels{i})));
            end
        end
        
        function p = defaultThresholds(p) %Need to update this with some better heuristic
            p.thresholds = zeros(0,numel(p.spotChannels));
            for i = 1:numel(p.spotChannels)
                p.thresholds(i) = max(1, round(mean(p.getIntensities(p.spotChannels{i}))));
            end
        end
        
        function p = makeCentroidList(p)
            p.centroidLists = cell(0, numel(p.spotChannels));
            for i = 1:numel(p.spotChannels)
                p.centroidLists{i} = sortrows(p.tabulateChannel(p.spotChannels{i}), 'GroupCount', 'descend');
                p.centroidLists{i}.expression_color = d2utils.expressionToColors(p.centroidLists{i}.GroupCount, p.expressionColorPal{p.paletteIdx});
                %Note that d2utils.expressionToColors assumed
                %min(GroupCount) = 0. If not, may need to update.
            end
        end
        
        function p = updateCentroidList(p, channel)
            channelIdx = ismember(p.spotChannels, channel);
            if any(channelIdx) % if the channel doesn't have channelType of FISH then it isn't in spotTable
                p.centroidLists{channelIdx}...
                    = sortrows(p.tabulateChannel(channel), 'GroupCount', 'descend');
                p.centroidLists{channelIdx}.expression_color = d2utils.expressionToColors(p.centroidLists{channelIdx}.GroupCount, p.expressionColorPal{p.paletteIdx});
            end
        end
        
        function p = updateExpressionColors(p)
            for i = 1:numel(p.spotChannels)
                p.centroidLists{i}.expression_color = d2utils.expressionToColors(p.centroidLists{i}.GroupCount, p.expressionColorPal{p.paletteIdx});
            end
        end
        
        function p = updateMaxDistance(p, distance)
            p.maxDistance = distance;
            p = updateAllSpotStatus(p);
        end
        
        function outTable = centroidTableInRect(p, channelIdx, rect)
%             channelIdx = ismember(p.spotChannels, channel);
            centroidIdx = p.centroidLists{channelIdx}.x >= rect(1)...
                & p.centroidLists{channelIdx}.x < rect(1) + rect(3) ...
                & p.centroidLists{channelIdx}.y >= rect(2)...
                & p.centroidLists{channelIdx}.y < rect(2) + rect(4);
            outTable = p.centroidLists{channelIdx}(centroidIdx,:);
        end
        
        function outTable = tabulateChannel(p, channel)
            idx = ismember(p.spots.channel, channel) & p.spots.status;
            tmpSpots = groupsummary(p.spots(idx, :), 'nearestNucID');
            outTable = outerjoin(p.nucleiObj.nuclei(p.nucleiObj.nuclei.status,{'nucID', 'x', 'y', 'colors'}), tmpSpots, 'Type', 'left', 'LeftKeys', 'nucID', 'RightKeys', 'nearestNucID');
            outTable{isnan(outTable.GroupCount), 'GroupCount'} = single(0);
        end
        
        function outTable = tabulateAllChannels(p)
            tmpSpots = groupsummary(p.spots(p.spots.status,:), {'channel', 'nearestNucID'}, 'IncludeEmptyGroups', true);
            spotLessCells = setdiff(p.nucleiObj.nuclei.nucID(p.nucleiObj.nuclei.status), tmpSpots.nearestNucID);
            nSpots = numel(p.spotChannels);
            spotLessTable = table(repelem(p.spotChannels', numel(spotLessCells), 1), repmat(spotLessCells, nSpots, 1), zeros(numel(spotLessCells)*nSpots, 1),...
                'VariableNames', tmpSpots.Properties.VariableNames);
            tmpSpots = [tmpSpots; spotLessTable];
            outTable = outerjoin(p.nucleiObj.nuclei(p.nucleiObj.nuclei.status,{'nucID', 'x', 'y'}), tmpSpots, 'Type', 'left', 'LeftKeys', 'nucID', 'RightKeys', 'nearestNucID');
            outTable.nearestNucID = [];
        end
        
        function updateScanSummary(p)
            scanSummaryTable = d2utils.parseScanSummary(p.scanObj.scanSummaryFile);
            scanSummaryTable{'spotChannels',:} = {strjoin(p.spotChannels)};
            scanSummaryTable{'thresholds',:} = {num2str(p.thresholds)};
            scanSummaryTable{'sigma',:} = {num2str(p.sigma)};
            scanSummaryTable{'aTrousMinThreshFactor',:} = {num2str(p.threshFactor)};
            scanSummaryTable{'minChannelIntensities',:} = {num2str(p.minChannelIntensities)};
            writetable(scanSummaryTable, p.scanObj.scanSummaryFile, 'WriteRowNames', true, 'WriteVariableNames', false, 'QuoteStrings', false, 'Delimiter', '\t')  
        end
        
        function p = parseScanSummary(p, inFileName)
            if isfile(inFileName)
                scanSummaryTable = d2utils.parseScanSummary(inFileName);
                if ismember('spotChannels', scanSummaryTable.Row)
                    p.spotChannels = split(scanSummaryTable{'spotChannels',1}{:})';
                    %Check if spot channels in the scan summary file match
                    %the channels in scanObject. If not, there may be
                    %issues later. Alternatively, don't read spotChannels
                    %from scan summary, just copy from scan object. 
                    if ~all(ismember(p.spotChannels, p.scanObj.channels(ismember(p.scanObj.channelTypes,{'FISH'}))))
                        error('Spot channels in scan summary file: %s do not match FISH channels in the scan object. You may need to edit or delete %s and try running again.', inFileName, inFileName)
                    end
                    fprintf('Setting %s as spot channels.\n', scanSummaryTable{'spotChannels',1}{:});
                else
                    p.spotChannels = p.scanObj.channels(ismember(p.scanObj.channelTypes,{'FISH'}));
                    fprintf('Setting %s as spot channels.\n', strjoin(p.spotChannels));
                end

                if ismember('thresholds', scanSummaryTable.Row)
                    p.thresholds = str2double(split(scanSummaryTable{'thresholds',1})');
                end
                
                if ismember('aTrousMinThreshFactor',scanSummaryTable.Row)
                    p.threshFactor= str2double(split(scanSummaryTable{'aTrousMinThreshFactor',1})');
                end
                
                if ismember('minChannelIntensities',scanSummaryTable.Row)
                    p.minChannelIntensities= str2double(split(scanSummaryTable{'minChannelIntensities',1})');
                end
                
                if ismember('sigma',scanSummaryTable.Row)
                    p.sigma= str2double(split(scanSummaryTable{'sigma',1})');
                end
                
            else
                fprintf('Unable to find %s in your current directory.\n', inFileName)
                p.spotChannels = p.scanObj.channels(ismember(p.scanObj.channelTypes,{'FISH'}));
                fprintf('Setting %s as spot channels.\n', strjoin(p.spotChannels));
            end
        end

        function saveSpotsTable(p)
            if ~isempty(p.spots)
                writetable(p.spots(:,{'spotID', 'x', 'y', 'intensity', 'nearestNucID', 'status', 'maskID', 'channel', 'distanceToNuc'}), p.spotsFile);
            else
                fprintf("spots table is empty. Run findSpots and try again")
            end
        end
        
        function exportSpotsSummary(p)
            if isempty(p.spots(p.spots.status,:)) || all(p.spots.nearestNucID == 0)
                fprintf("There are no valid spots or spots are not assigned to cells. Try running findSpots and assignSpotsToNuclei")
            else
                outTable = p.tabulateAllChannels;
                writetable(outTable, p.spotsSummaryFile);
            end
        end
        

        function userInputThresholds(p,userInputThresholds)
            if numel(userInputThresholds)~=numel(p.spotChannels)
                error('user-input thresholds is 1x%i but should be 1x%i, corresponding to %i spot channels',numel(userInputThresholds),numel(p.spotChannels),numel(p.spotChannels))
            end
            p.thresholds=userInputThresholds;
        end
        

    end
    
end 