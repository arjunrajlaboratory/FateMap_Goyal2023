classdef d2ThresholdController < handle
    properties (Access = public)
        
        %Model objects
        spotTable
        scanObj
        maskObj
        nucleiObj
        scanDim
        viewObj
        threshZoom
        viewRect %Image coords not axes coords
        panRect
        cellViewRadius = 750
        zoomROI
        zoomStart
        zoomRect
        fixedZoom = false
        panStart
        channelIdx
        zoomMode = true
        scatterH %Not sure if we need these handle. 
        imageH
        spotScatterH
        nucleiScatterH
        maskH
        handPointStruct 
        imagesInView
        dapiInView
        
        resizedImg %Not sure if we want to save these or regenerate
        resizedDapi
        
    end
    
    methods
        function p = d2ThresholdController(view, scanObj, spotTable, maskObj, nucleiObj)
            p.viewObj = view;
            p.spotTable = spotTable;
            p.scanObj = scanObj;
            p.maskObj = maskObj;
            p.nucleiObj = nucleiObj;
            p.scanDim = size(p.scanObj.dapiStitch);
            p.viewRect = [1, 1, p.scanDim];
            p.startup()
        end
        
        function startup(p)
            p.channelIdx = p.viewObj.channelPopup.Value;
            p.spotTable.makeIntensitiesToPlot();
            p.plotIntensityHistogram()
            p.plotIntensityThreshold()
%             p.updateColorMap();  %Specifying colors directly in scatter instead
            p.plotScatterMain();
            p.updateImageInView();
            p.handPointStruct = d2utils.makePointStruct();
            %p.updateMainAxes()
        end
        
        function p = changeChannel(p, ~, ~)
            p.channelIdx = p.viewObj.channelPopup.Value;
            %Update centroid listbox
            p.updateCentroidListView();
%             p.updateColorMap(); %Specifying colors directly in scatter instead
            p.updateMainAxes();
            p.plotIntensityHistogram();
        end
        
        function p = updateMainAxes(p, ~, ~)
            if isvalid(p.viewObj.mainAxes.Children)
                delete(get(p.viewObj.mainAxes, 'Children'));
            end
            
            if logical(p.viewObj.scatterCheckBox.Value)
                p.plotScatterMain();
            else
                p.showImage();
                p.overlaySpots();
                p.overlayNuclei();
                p.overlayMasks();
            end
        end
        
        function updateCentroidListView(p)
            p.viewObj.centroidList.String = string(p.spotTable.centroidLists{p.channelIdx}.GroupCount);
        end
        
        function p = centroidSelected(p, ~, ~)
            if strcmp(get(p.viewObj.figHandle, 'SelectionType'), 'open') %Respond to double mouse click
                cellIdx = get(p.viewObj.centroidList, 'Value');
                cellPos = p.spotTable.centroidLists{p.channelIdx}{cellIdx, {'x', 'y'}};
                p.viewRect = p.getRectAroundPoint(cellPos);
                p.updateImageInView();
                p.updateMainAxes();
            end
        end
        
        function outRect = getRectAroundPoint(p, point)
            startPos = max([1, 1], point - p.cellViewRadius); %Avoid rect outside range of scan. 
            startPos = min(size(p.scanObj.dapiStitch)- p.cellViewRadius+1,startPos);
            outRect = [startPos, [2, 2] * p.cellViewRadius];
        end
        
        function plotIntensityHistogram(p)
            %channelIdx = p.viewObj.channelPopup.Value;
            if ~isempty(p.viewObj.histogramLineH)
                delete(p.viewObj.histogramLineH)
            end
            intensities = p.spotTable.intensitiesToPlot{p.channelIdx};
            logRank = log(numel(intensities):-1:1);
            p.viewObj.histogramLineH = line(p.viewObj.threshAxes, intensities, logRank, ...
                'HitTest', 'off', ...
                'Color', 'k');
            
            xAxisMin = intensities(1);
            xAxisMax = intensities(end) * 1.05;
            
            set(p.viewObj.threshAxes, 'XLim', [xAxisMin xAxisMax]);
            
            yaxismax = logRank(1)*1.1;
          
            set(p.viewObj.threshAxes, 'YLim', [0 yaxismax]);
        end
        
        function plotIntensityThreshold(p)
            %channelIdx = p.viewObj.channelPopup.Value;
            if ~isempty(p.viewObj.thresholdLineH) && ishandle(p.viewObj.thresholdLineH)
                delete(p.viewObj.thresholdLineH)
            end
            if isempty(p.spotTable.thresholds)
                p.spotTable.defaultThresholds();
            end
            threshold = p.spotTable.thresholds{p.channelIdx};
            yaxis = get(p.viewObj.threshAxes, 'Ylim');
            p.viewObj.thresholdLineH = line(p.viewObj.threshAxes, [threshold threshold], yaxis,...
                'Color', 'b', 'HitTest', 'off');
            p.viewObj.threshValue.String = num2str(threshold);
        end
        
        function plotScatterMain(p)
            %Could probably save on time by adding colormap to
            %centroidTable in the spotTableObject.
            centroidsInView = p.spotTable.centroidTableInRect(p.channelIdx, p.viewRect);
            centroidsInView = flipud(centroidsInView); %In order to plot cells with 0 expression in the back
            set(p.viewObj.mainAxes, 'XLim', [p.viewRect(2),  p.viewRect(2)+p.viewRect(4)])
            set(p.viewObj.mainAxes, 'YLim', [p.viewRect(1),  p.viewRect(1)+p.viewRect(3)])
            hold(p.viewObj.mainAxes, 'on')
            %Note that it is slightly slower specifying RGB colors with expression_color
            %instead of changing the colormap on the axes however, specifying RGB colors
            %ensures that the colors don't change with the view 
            p.scatterH = scatter(centroidsInView.y, centroidsInView.x,...
                20, centroidsInView.expression_color, 'filled',...
                'Parent', p.viewObj.mainAxes, 'HitTest','off');
            hold(p.viewObj.mainAxes, 'off')
            %colorbar(p.viewObj.mainAxes, 'Location', 'eastoutside','HitTest','off')
        end
        
        function zoomInPressed(p, ~, ~)
            p.zoomMode = true;
            %iptPointerManager(p.viewObj.figHandle, 'disable');
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
            %Below are old functions that draw rectangle using matlab
            %ROI instead of mainAxesZoom fcn. Probably slightly more efficient but doesn't allow using
            %the console until roi completed. 
            %set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDownZoomMode}) %need a zoomMode that doesn't use 'normal' clicks
%             while p.zoomMode
%                 if p.fixedZoom
%                     p.zoomROI = drawrectangle(p.viewObj.mainAxes, 'Color', 'r', 'FaceAlpha', 0, 'FixedAspectRatio', true);
%                 else 
%                     p.zoomROI = drawrectangle(p.viewObj.mainAxes, 'Color', 'r', 'FaceAlpha', 0);
%                 end
%                 
%                 if isempty(p.zoomROI.Position) %if escape pressed
%                     break
%                 elseif all(p.zoomROI.Position(3:4) > 4) %minimum view size
%                     p.viewRect  = d2utils.coordToPixelRect(round(p.zoomROI.Position));
%                     delete(p.zoomROI) %Possibly unnecessary since Children deleted in p.updateMainAxes(); 
%                     p.updateImageInView();
%                     p.updateMainAxes();
%                     p.fixedZoom = false;
%                 end
%             end
        end
    
        function figWindowDown(p, ~, ~)
            if d2utils.pointInSideRect(p.viewRect, get(p.viewObj.mainAxes, 'CurrentPoint'))
                switch(p.getSelectionType)
                    %'normal' click was used for drawing rectangle before using
                    %rectangle roi
                    case 'normal'
                        if p.zoomMode
                            p.zoomStart = get(p.viewObj.mainAxes, 'CurrentPoint');
                            set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopDragFcn});
                            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.mainAxesZoom});
                        else
                            p.panStart =  get(p.viewObj.mainAxes, 'CurrentPoint');
                            set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopPan})
                            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.panView});
                        end
                    case 'open'
                        p.viewRect = [1, 1, p.scanDim];
                        p.updateImageInView();
                        p.updateMainAxes();
                    case 'alt'
                        p.viewRect = d2utils.expandView2x(p.viewRect, p.scanDim);
                        p.updateImageInView();
                        p.updateMainAxes();
                end
            end
        end
        
        function mainAxesZoom(p, ~, ~)
            if ~isempty(p.viewObj.zoomH)
                delete(p.viewObj.zoomH)
            end
            
            if p.fixedZoom
                p.getSelectedRectangleCoordsFixed()
            else
                p.getSelectedRectangleCoords()
            end
            
            if all(p.zoomRect(3:4) > 4) %min view size
                p.viewObj.zoomH = rectangle(...
                    'Position', p.zoomRect, ...
                    'EdgeColor', 'r', ...
                    'Parent', p.viewObj.mainAxes,...
                    'Hittest', 'off');
            end

        end
        
        function stopDragFcn(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            if ~isempty(p.zoomRect) && all(p.zoomRect > 0)
                p.viewRect  = d2utils.coordToPixelRect(p.zoomRect); % Should notify event "viewChange"
                %Update imagesInView
                p.updateImageInView();
                %Update main axes
                p.updateMainAxes();
                p.zoomRect = [];
            end
            delete(p.viewObj.zoomH)
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '');
        end
        
        function getSelectedRectangleCoords(p)
            currentPoint = get(p.viewObj.mainAxes, 'CurrentPoint');
            p.zoomRect = round([min(p.zoomStart(1,1:2), currentPoint(1,1:2)),...
                abs(p.zoomStart(1,1:2) - currentPoint(1,1:2))]);
        end
        
        function getSelectedRectangleCoordsFixed(p)
            currentPoint = get(p.viewObj.mainAxes, 'CurrentPoint');
            sz = max(abs(p.zoomStart(1,1:2) - currentPoint(1,1:2)));
            p.zoomRect = round([min(p.zoomStart(1,1:2), currentPoint(1,1:2)), [sz sz]]);
        end
        
        function panViewPressed(p, ~, ~)
            if any(p.viewRect(3:4) < p.scanDim) %No need to pan if view is of entire scan
                p.zoomMode = false; 
                %iptPointerManager(p.viewObj.figHandle, 'enable');
                %iptSetPointerBehavior(p.viewObj.mainAxes, @(hfig, currentPoint) set(hfig, 'Pointer', 'hand'));
                set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
            else
                p.zoomMode = false; 
                set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            end
        end
                
        function stopPan(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '')
        end
        
        function panView(p, ~, ~)
            currentPoint = get(p.viewObj.mainAxes, 'CurrentPoint');
            displacement = p.panStart(1,1:2) - currentPoint(1,1:2);
            p.viewRect = d2utils.updateViewPanning(p.viewRect, displacement, p.scanDim);
            p.updateImageInView;
            p.updateMainAxes;
        end
        
        function thresholdButtonDown(p, ~, ~)
            %set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopThresholdDrag});
            %set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.thresholdDrag});
            currentPoint = get(p.viewObj.threshAxes, 'CurrentPoint');
            if abs(currentPoint(1,1) - p.viewObj.thresholdLineH.XData(1)) < 150
                set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopThreshDrag});
                set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.dragThresh});
            end
        end
        
        function dragThresh(p, ~, ~)
            currentPont = get(p.viewObj.threshAxes, 'CurrentPoint');
            newThresh = max(currentPont(1,1), 1);
            newThresh = min(newThresh, p.viewObj.threshAxes.XLim(2));
            if isvalid(p.viewObj.thresholdLineH)
                set(p.viewObj.thresholdLineH, 'XData', [newThresh, newThresh])
                p.viewObj.threshValue.String = num2str(round(newThresh));
            else
                yaxis = get(p.viewObj.threshAxes, 'Ylim');
                p.viewObj.thresholdLineH = line(p.viewObj.threshAxes, [newThresh, newThresh], yaxis,...
                    'Color', 'b', 'HitTest', 'off');
                p.viewObj.threshValue.String = num2str(round(newThresh));
            end
        end
        
        function stopThreshDrag(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            %Should notify threshold change
            newThresh = str2double(p.viewObj.threshValue.String);
            p.spotTable.setThreshold(p.spotTable.spotChannels{p.channelIdx}, newThresh);
            p.viewObj.centroidList.String = string(p.spotTable.centroidLists{p.channelIdx}.GroupCount);
%             p.updateColorMap(); %Specifying colors directly in scatter instead
            p.updateMainAxes();
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '');
        end
        
        function updateColorMap(p)
            %Can use this function to change the axes colormap instead of
            %specifying colors in the scatter function. However, if relying
            %on the axes colormap, then the point colors may change
            %depending on the field of view (specifically the expression level of other cells in the FOV).
            spotCounts = p.spotTable.centroidLists{p.channelIdx};
            colors = d2utils.expressionToColors(min(spotCounts.GroupCount):max(spotCounts.GroupCount));
            set(p.viewObj.mainAxes, 'Colormap', colors);
        end
        
        function threshValueChange(p, ~, ~)
            if ~isnan(str2double(p.viewObj.threshValue.String)) && isreal(str2double(p.viewObj.threshValue.String))
                newThresh = round(str2double(p.viewObj.threshValue.String));
                set(p.viewObj.thresholdLineH, 'XData', [newThresh, newThresh])
                %Update threshold for spots
                p.spotTable.setThreshold(p.spotTable.spotChannels{p.channelIdx}, newThresh);
                p.viewObj.centroidList.String = string(p.spotTable.centroidLists{p.channelIdx}.GroupCount);
%                 p.updateColorMap(); %Specifying colors directly in scatter instead
                p.updateMainAxes();
            else
                %If not number, return to previous threshold or default
                %threshold
                if isvalid(p.viewObj.thresholdLineH)
                    oldThresh = get(p.viewObj.thresholdLineH, 'XData');
                    set(p.viewObj.threshValue, 'String', num2str(oldThresh(1)))
                else
                    threshold = p.spotTable.thresholds{p.channelIdx};
                    yaxis = get(p.viewObj.threshAxes, 'Ylim');
                    p.viewObj.thresholdLineH = line(p.viewObj.threshAxes, [threshold threshold], yaxis,...
                        'Color', 'b', 'HitTest', 'off');
                    p.viewObj.threshValue.String = num2str(threshold);
                    p.spotTable.setThreshold(p.spotTable.spotChannels{p.channelIdx}, num2str(threshold));
%                     p.updateColorMap(); %Specifying colors directly in scatter instead
                    p.updateMainAxes();
                end
            end
            
        end
%         function threshValueKey(p, ~, evt) %Could use this to toggle up and down with arrow keys
%             if strcmp(evt.Key, 'return')
%                 disp(p.viewObj.threshValue.String)
%             end
%         end
       
        function p = updateImageInView(p)
            p.imagesInView = cell(0, numel(p.spotTable.spotChannels));
            if p.viewRect(3) * p.viewRect(4) < 4000001
                for i = 1:numel(p.spotTable.spotChannels)
                    p.imagesInView{i} = p.scanObj.getImageRect(p.spotTable.spotChannels{i}, p.viewRect);
                end
                p.dapiInView = p.scanObj.getDapiImage(p.viewRect);
            elseif p.viewRect(3) * p.viewRect(4) < 64000001
                for i = 1:numel(p.spotTable.spotChannels)
                    p.imagesInView{i} = p.scanObj.getSmallImageRect(p.spotTable.spotChannels{i}, p.viewRect);
                end
                p.dapiInView = p.scanObj.getSmallDapiImage(p.viewRect);
            else
                disp('View is too large to display image. Plotting scatter')
                set(p.viewObj.scatterCheckBox, 'Value', 1)
            end
        end
        
        function showImage(p)
           %Adjust contrast
            contrastIn = [get(p.viewObj.lowerContrastSlider, 'Value'), get(p.viewObj.upperContrastSlider, 'Value')];
            tmpIm = imadjust(p.imagesInView{p.channelIdx}, contrastIn, []);
            %Decide if overlay with DAPI?
            tmpRGB = cat(3, tmpIm, tmpIm, tmpIm+p.dapiInView);
            xlimits = [p.viewRect(2),  p.viewRect(2)+p.viewRect(4)];
            ylimits = [p.viewRect(1),  p.viewRect(1)+p.viewRect(3)];
            %axis(p.viewObj.mainAxes, [xlimits, ylimits], 'square')
            set(p.viewObj.mainAxes, 'XLim', xlimits)
            set(p.viewObj.mainAxes, 'YLim', ylimits)
            hold(p.viewObj.mainAxes, 'on')
            p.imageH = imshow(tmpRGB, 'XData', xlimits, 'YData', ylimits, 'Parent', p.viewObj.mainAxes);
            %axis fill
            %pbaspect auto
            set(p.viewObj.mainAxes, 'Visible', 'on')
            hold(p.viewObj.mainAxes, 'off')
        end
        
        function overlaySpots(p, ~, ~)
            if logical(p.viewObj.spotsCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                [spotsInView, spotIdx] = p.spotTable.getValidSpotsInRect(p.spotTable.spotChannels{p.channelIdx}, p.viewRect);
                hold(p.viewObj.mainAxes, 'on')
                p.spotScatterH = scatter(spotsInView.y, spotsInView.x, 10, spotsInView.colors,...
                    'Parent', p.viewObj.mainAxes, 'HitTest','off');
                hold(p.viewObj.mainAxes, 'off')
            else
                delete(p.spotScatterH)
            end
        end
        
        function overlayNuclei(p, ~, ~)
            if logical(p.viewObj.centroidsCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                outTableTmp = p.spotTable.centroidTableInRect(p.channelIdx, p.viewRect);
                hold(p.viewObj.mainAxes, 'on')
                p.nucleiScatterH = scatter(outTableTmp.y, outTableTmp.x, 30, outTableTmp.colors, 'filled',...
                    'Parent', p.viewObj.mainAxes, 'HitTest','off');
                hold(p.viewObj.mainAxes, 'off')
            else
                delete(p.nucleiScatterH)
            end
        end
        
        function overlayMasks(p, ~, ~)
            masks = findobj(p.viewObj.mainAxes, 'type', 'images.roi.freehand');
            delete(masks);
            if logical(p.viewObj.masksCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                maskTableTmp = p.maskObj.getChannelMasksInRect(p.viewRect, p.spotTable.spotChannels{p.channelIdx});
                maskIDs = unique(maskTableTmp.maskID);
                maskIDs(maskIDs == 0) = [];
                for i = 1:numel(maskIDs)
                    drawfreehand(p.viewObj.mainAxes, 'Position', maskTableTmp{maskTableTmp.maskID == maskIDs(i), {'y', 'x'}},...
                        'Color', 'red', 'InteractionsAllowed', 'none');
                end
                cellMasksTmp = p.maskObj.getChannelMasksInRect(p.viewRect, 'dapi');
                maskIDs = unique(cellMasksTmp.maskID);
                maskIDs(maskIDs == 0) = [];
                for i = 1:numel(maskIDs)
                    drawfreehand(p.viewObj.mainAxes, 'Position', cellMasksTmp{cellMasksTmp.maskID == maskIDs(i), {'y', 'x'}},...
                        'Color', 'blue', 'InteractionsAllowed', 'none');
                end
            end
        end
        
        function scatterCallback(p, ~, ~)
            if p.viewObj.scatterCheckBox.Value == 0
                p.updateImageInView();
            end
            p.updateMainAxes();
        end
        
        function addSpotMask(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'off')
            channel = p.spotTable.spotChannels{p.channelIdx};
            p.maskH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            %addlistener(p.maskH, 'DrawingFinished', @p.maskSpots);
            if ~isempty(p.maskH.Position) %Allows 'escape' from ROI
                p.maskSpots(p.maskH, channel)
            end
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function maskSpots(p, roi, channel)
            %tmpPoly = roi.Position;
            set(p.viewObj.masksCheckBox, 'Value', true)
            p.maskObj.addMaskLocalCoords(roi.Position, channel);
            delete(roi)
            tic
            p.spotTable.addNewMask(channel);
            p.spotTable.updateSpotStatus(channel);
            toc
            p.spotTable.updateCentroidList(channel);
            p.updateCentroidListView();
            p.updateMainAxes();
        end
        
        function addCellMask(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'off')
            p.maskH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            if ~isempty(p.maskH.Position) %Allows 'escape' from ROI
                p.maskCells(p.maskH)
            end
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function maskCells(p, roi)
            %tmpPoly = roi.Position;
            channel = p.spotTable.spotChannels{p.channelIdx};
            set(p.viewObj.masksCheckBox, 'Value', true)
            p.maskObj.addMaskLocalCoords(roi.Position, 'dapi');
            delete(roi)
            p.nucleiObj.addNewMask();
            p.spotTable.updateSpotStatus(channel);
            p.spotTable.updateCentroidList(channel);
            p.updateCentroidListView();
            p.updateMainAxes();
        end
        
        function deleteMask(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'off')
            channel = p.spotTable.spotChannels{p.channelIdx};
            [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.   
            if ~isempty(x)
                ptsInView = d2utils.getPtsInsideView([x, y], p.viewRect);
                p.maskObj.removeMasksByLocalPoints(ptsInView, p.viewRect);
                tic
                p.nucleiObj.removeMasks();
                toc
                %p.nucleiObj.updateMasksInRect(p.viewRect);
                tic
                p.spotTable.removeMasks2(channel, p.viewRect);
                toc
                %p.spotTable.updateMasksInRect(channel, p.viewRect);
                p.spotTable.updateCentroidList(channel);
                p.updateCentroidListView();
                p.updateMainAxes();
                set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'on')
            end
        end
        
        function addCells(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'off')
            [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.   
            if ~isempty(x)
                ptsInView = d2utils.getPtsInsideView([x, y], p.viewRect);
                p.nucleiObj.addCell(ptsInView(:,2), ptsInView(:,1));
                p.spotTable.assignSpotsInRect(p.viewRect)
                %p.nucleiObj.updateMasksInRect(p.viewRect);
                tic
                p.spotTable.updateAllSpotStatus();
                toc
                tic
                p.spotTable.makeCentroidList();
                toc
                p.updateCentroidListView();
                p.updateMainAxes();
                set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'on')
            end
        end
        
        function deleteCells(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes], 'Enable', 'off')
            [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.
            if ~isempty(x)
                ptsInView = d2utils.getPtsInsideView([x, y], p.viewRect);
            end
        end
        
        function type = getSelectionType(p)
            type = get(p.viewObj.figHandle, 'SelectionType');
        end
                
        function keyPressFcns(p, ~, evt)
            keyPressed = evt.Key;
            switch(keyPressed)
                case 'z'
                    p.zoomInPressed();
                case 'p'
                    p.panViewPressed();
                case 'x'
                    set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
                    %iptPointerManager(p.viewObj.figHandle, 'disable');
                case 'shift'
                    disp('shift')
                    p.fixedZoom = ~p.fixedZoom;
%                 case 'return' %'return' key is used to end mask polygon.
%                 % If want to use 'return' here, then change mask callback.  
%                 
%                     if isvalid(p.zoomROI)
%                         p.viewRect = d2utils.coordToPixelRect(p.zoomROI.Position);
%                         delete(p.zoomROI)
%                         delete(p.scatterH)
%                         set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown});
%                         p.updateImageInView
%                         Update main axes
%                         p.updateMainAxes();
%                         drawrectangle('Position', pos, 'Color', 'r', 'InteractionsAllowed', 'none');
%                     end
                    
                
                    
            end
                
        end
        
        
        
        
        

        
        
    end
    
end