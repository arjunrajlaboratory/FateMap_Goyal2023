classdef d2IFController < handle
    properties (Access = public)
        
        %Model objects
        scanObj
%         maskObj
        IFboundaries
        IFtable
        
        %View
        viewObj
        
        thumbCntrlr
        
        viewRect %Image coords not axes coords
        panRect
        cellViewRadius = 500
        zoomROI
        zoomStart
        zoomRect
        quantMetricDict
        quantMetric
        fixedZoom = false
        panStart
        channelIdx
        zoomMode = true
        deleteROI = false
        scatterH %Not sure if we need these handles. 
        imageH
        nucleiPlotH
        cellPlotH
        holesMaskH
        imgMaskH
        cellMaskH
        nucH
        cellH
        maskH
        imagesInView
        dapiInView
        
        numCells
    end
    
    methods
        function p = d2IFController(view, scanObj, IFboundaries, IFtable)
            p.viewObj = view;
            p.scanObj = scanObj;
%             p.maskObj = maskObj;
            p.IFboundaries = IFboundaries;
            p.IFtable = IFtable;
            p.viewRect = [1, 1, min(p.scanObj.stitchDim, [25000 25000])];
            p.startup()
        end
        
        function startup(p)
            p.channelIdx = p.viewObj.channelPopup.Value;
            p.updateImageInView();
            p.quantMetricDict = containers.Map({'mean nucleus', 'mean cytoplasm', 'sum nucleus', 'sum cytoplasm'},...
                            {'meanNuc', 'meanCyto', 'sumNuc', 'sumCyto'}); %Hard-coded for now. Can make this more functional in the future
            p.quantMetric = get(p.viewObj.quantMetric.SelectedObject, 'String');
%             p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
            p.numCells = size(p.IFtable.centroidLists{p.channelIdx}, 1);
            if isempty(p.IFtable.centroidLists{p.channelIdx})
                p.showImage();
                p.overlayMasks();
            else
                p.plotScatterMain();
            end
        end
        
        function p = changeChannel(p, ~, ~)
            p.channelIdx = p.viewObj.channelPopup.Value;
            set(p.viewObj.IFCheckBox, 'String', sprintf('%s (f)', p.IFtable.channels{1}));
            %Update centroid listbox
            p.updateCentroidListView();
            p.updateMainAxes();
        end
        
        function changeColormap(p, ~, ~)
            p.IFtable.paletteIdx = p.viewObj.colormapPopup.Value;
            p.IFtable.updateExpressionColors(p.quantMetricDict(p.quantMetric)); 
            p.updateMainAxes();
        end
        
        function p = updateMainAxes(p, ~, ~)
            if isvalid(p.viewObj.mainAxes.Children)
                delete(get(p.viewObj.mainAxes, 'Children'));
            end
            
            if logical(p.viewObj.scatterCheckBox.Value)
                p.plotScatterMain();
            else
                p.showImage();
                p.overlayNuclei();
                p.overlayCells();
                p.overlayMasks();
            end
        end
        
        function updateCentroidListView(p)
            set(p.viewObj.centroidList, 'String', string(p.IFtable.centroidLists{p.channelIdx}{:,4})); %  %Add var (e.g. meanNuc, meanCyto etc.)
            p.numCells = size(p.IFtable.centroidLists{p.channelIdx}, 1); %Although this value doesn't depend on the channel. Easy to put this here rather than everywhere where cell # can change. 
        end
        
        function p = centroidSelected(p, ~, ~)
            if strcmp(get(p.viewObj.figHandle, 'SelectionType'), 'open') %Respond to double mouse click
                cellIdx = get(p.viewObj.centroidList, 'Value');
                cellPos = p.IFtable.centroidLists{p.channelIdx}{cellIdx, {'x', 'y'}};
                p.viewRect = d2utils.getRectAroundPoint(cellPos, 2 * p.cellViewRadius, 2 * p.cellViewRadius, p.scanObj.stitchDim);
                set(p.viewObj.scatterCheckBox, 'Value', 0)
                p.updateImageInView();
                p.updateMainAxes();
                p.thumbCntrlr.overlayThumbnailRect();
            end
        end

        function quantMetricChanged(p, ~, evt)
            p.quantMetric = evt.NewValue.String;
            cellIdx = get(p.viewObj.centroidList, 'Value');
            cellID = p.IFtable.centroidLists{p.channelIdx}.cellID(cellIdx);
            p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
            p.updateCentroidListView();
            newCellIdx = find(p.IFtable.centroidLists{p.channelIdx}.cellID == cellID);
            set(p.viewObj.centroidList, 'Value', newCellIdx);
            p.updateMainAxes();
        end
        
        function selectionToolChanged(p, ~, evt)
            switch(evt.NewValue.String)
                case 'draw'
                    p.deleteROI = true;
                case 'points' 
                    p.deleteROI = false;
            end
        end
        
        function plotScatterMain(p)
            centroidsInView = p.IFtable.centroidTableInRect(p.channelIdx, p.viewRect);
            centroidsInView = flipud(centroidsInView); %In order to plot cells with 0 expression in the back
            set(p.viewObj.mainAxes, 'XLim', [p.viewRect(2),  p.viewRect(2)+p.viewRect(4)])
            set(p.viewObj.mainAxes, 'YLim', [p.viewRect(1),  p.viewRect(1)+p.viewRect(3)])
            hold(p.viewObj.mainAxes, 'on')
            %Note that it is slightly slower to specify RGB colors with expression_color
            %instead of changing the colormap on the axes however, specifying RGB colors
            %ensures that the colors don't change with the view 
            p.scatterH = scatter(centroidsInView.y, centroidsInView.x,...
                30, centroidsInView.expression_color, 'filled',...
                'Parent', p.viewObj.mainAxes, 'HitTest','off');
            hold(p.viewObj.mainAxes, 'off')
            %colorbar(p.viewObj.mainAxes, 'Location', 'eastoutside','HitTest','off')
        end
        
        function zoomInPressed(p, ~, ~)
            p.zoomMode = true;
            %iptPointerManager(p.viewObj.figHandle, 'disable');
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
    
        function figWindowDown(p, ~, ~)
            if d2utils.pointInSideViewRect(p.viewRect, get(p.viewObj.mainAxes, 'CurrentPoint'))
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
                        p.viewRect = [1, 1, min(p.scanObj.stitchDim, [25000 25000])];
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.thumbCntrlr.overlayThumbnailRect();
                    case 'alt'
                        p.viewRect = d2utils.expandView2x(p.viewRect, p.scanObj.stitchDim);
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.thumbCntrlr.overlayThumbnailRect();
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
                p.thumbCntrlr.overlayThumbnailRect();
                p.zoomRect = [];
                p.fixedZoom = false;
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
            if any(p.viewRect(3:4) < p.scanObj.stitchDim) %No need to pan if view is of entire scan
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
            p.viewRect = d2utils.updateViewPanning(p.viewRect, displacement, p.scanObj.stitchDim);
            p.updateImageInView;
            p.updateMainAxes;
            p.thumbCntrlr.overlayThumbnailRect();
        end
        
        function p = updateImageInView(p)
            p.imagesInView = cell(0, numel(p.IFtable.channels));
            if prod(p.viewRect(3:4)) < 4000001
                for i = 1:numel(p.IFtable.channels)
                    p.imagesInView{i} = p.scanObj.getImageRect(p.IFtable.channels{i}, p.viewRect);
                end
                p.dapiInView = p.scanObj.getDapiImage(p.viewRect);
            elseif prod(p.viewRect(3:4)) < 64000001
                for i = 1:numel(p.IFtable.channels)
                    p.imagesInView{i} = p.scanObj.getSmallImageRect(p.IFtable.channels{i}, p.viewRect);
                end
                p.dapiInView = p.scanObj.getSmallDapiImage(p.viewRect);
            else
                disp('View is too large to display image. Plotting scatter')
                set(p.viewObj.scatterCheckBox, 'Value', 1)
            end
        end
        
        function showImage(p)
            if get(p.viewObj.IFCheckBox, 'Value') 
                contrastIn = [get(p.viewObj.lowerContrastSlider, 'Value'), get(p.viewObj.upperContrastSlider, 'Value')];
                tmpIm = imadjust(p.imagesInView{p.channelIdx}, contrastIn, []);
%                 tmpRGB = cat(3, tmpIm, tmpIm, tmpIm+p.dapiInView);
            else
                tmpIm = zeros(size(p.imagesInView{p.channelIdx}), 'uint16');
            end
            
            if get(p.viewObj.dapiCheckBox, 'Value')
                tmpRGB = cat(3, tmpIm, tmpIm, tmpIm+p.dapiInView);
            else
                tmpRGB = repelem(tmpIm, 1, 1, 3);
            end
            
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
        
        function overlayCells(p, ~, ~)
            %Unlike overlayNuclei, overlayCells does not plot disconnected
            %boundaries. To plot disconnected cells, change
            %polyvect2patchMats to polyvect2patchMats2 as in overlayNuclei.
            if logical(p.viewObj.cellBordersCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                cellsInView = p.IFboundaries.getCellBoundariesInRect(p.IFboundaries.channels(p.channelIdx), p.viewRect);
                [cellXmat, cellYmat] = d2utils.polyvect2patchMats(cellsInView.cellBoundary); 
                hold(p.viewObj.mainAxes, 'on')
                p.cellPlotH = patch(cellYmat, cellXmat, reshape(cellsInView.colors, [], 1, 3),...
                    'FaceAlpha', 0.3, 'Parent', p.viewObj.mainAxes, 'HitTest','off', 'Tag', 'NucBoundary');
%                 for i = 1:size(cellsInView, 1)
%                     patch(cellsInView.cellBoundary(i).Vertices(:,2), cellsInView.cellBoundary(i).Vertices(:,1), cellsInView.colors(i,:),...
%                     'FaceAlpha', 0.3, 'Parent', p.viewObj.mainAxes, 'HitTest','off', 'Tag', 'CellBoundary');
%                 end
                hold(p.viewObj.mainAxes, 'off')
            else
                delete(p.cellPlotH)
%                 delete(findobj(p.viewObj.mainAxes, 'Tag', 'CellBoundary')) 
            end
        end
        
        function overlayNuclei(p, ~, ~)
            if logical(p.viewObj.nucleiBordersCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                nucleiInView = p.IFboundaries.getNucBoundariesInRect(p.IFboundaries.channels(p.channelIdx), p.viewRect);
%                 tic
%                 noHolesIdx = [nucleiInView.nucBoundary.NumHoles] == 0; Old strategy for plotting nuclei with holes
%                 [nucXmat, nucYmat, colorMat] = d2utils.polyvect2patchMats2(nucleiInView.nucBoundary(noHolesIdx), nucleiInView.colors(noHolesIdx, :));
                [nucXmat, nucYmat, colorMat] = d2utils.polyvect2patchMats2(nucleiInView.nucBoundary, nucleiInView.colors);
                hold(p.viewObj.mainAxes, 'on')
                p.nucleiPlotH = patch(nucYmat, nucXmat, reshape(colorMat, [], 1, 3),...
                    'FaceAlpha', 0.3, 'Parent', p.viewObj.mainAxes, 'HitTest','off', 'Tag', 'NucBoundary');
%                 if any(~noHolesIdx) Old strategy for plotting nuclei with holes
%                     tmpNuc = nucleiInView(~noHolesIdx, :)
%                     for i = 1:size(tmpNuc, 1)
%                         tmpPoly = polyshape(fliplr(tmpNuc.nucBoundary(i).Vertices));
%                         plot(tmpPoly, 'FaceAlpha', 0.3, 'FaceColor', tmpNuc.colors(i,:), 'Parent', p.viewObj.mainAxes, 'HitTest','off', 'Tag', 'NucBoundary');
%                     end
%                 end
%                 holesTmp = p.IFtable.maskObj.getChannelMasksInRect(p.viewRect, 'dapi');
%                 if ~isempty(holesTmp)
%                     [holesXmat, holesYmat] = d2utils.masktable2patchMats(holesTmp);
%                     hold(p.viewObj.mainAxes, 'on')
%                     p.holesMaskH = patch(holesYmat, holesXmat, 'white',...
%                         'FaceAlpha', 0.7, 'Parent', p.viewObj.mainAxes, 'HitTest','off', 'Tag', 'holesMask');
%                     hold(p.viewObj.mainAxes, 'off')
%                 end
                hold(p.viewObj.mainAxes, 'off')
%                 toc
            else
                delete(p.nucleiPlotH)
%                 delete(findobj(p.viewObj.mainAxes, 'Tag', 'NucBoundary'))
            end
        end
        
        function overlayMasks(p, ~, ~)
            if logical(p.viewObj.masksCheckBox.Value) && ~logical(p.viewObj.scatterCheckBox.Value)
                imgMasksTmp = p.IFtable.maskObj.getChannelMasksInRect(p.viewRect, p.IFtable.channels{p.channelIdx});
%                 imgMasksTmp(imgMasksTmp.maskID == 0,:) = [];
                if ~isempty(imgMasksTmp)
                    [imgMaskXmat, imgMaskYmat] = d2utils.masktable2patchMats(imgMasksTmp);
                    hold(p.viewObj.mainAxes, 'on')
                    p.imgMaskH = patch(imgMaskYmat, imgMaskXmat, 'white',...
                        'FaceAlpha', 0.7, 'Parent', p.viewObj.mainAxes, 'HitTest','off', 'Tag', 'imgMask');
                    hold(p.viewObj.mainAxes, 'off')
                end
                cellMasksTmp = p.IFboundaries.maskObj.getChannelMasksInRect(p.viewRect, p.IFboundaries.channels{p.channelIdx});
%                 cellMasksTmp(cellMasksTmp.maskID == 0,:) = [];
                if ~isempty(cellMasksTmp)
                    [cellMaskXmat, cellMaskYmat] = d2utils.masktable2patchMats(cellMasksTmp);
                    hold(p.viewObj.mainAxes, 'on')
                    p.cellMaskH = patch(cellMaskYmat, cellMaskXmat, [0.6 0.6 0.6],...
                        'FaceAlpha', 0.7, 'Parent', p.viewObj.mainAxes, 'HitTest','off', 'Tag', 'cellMask');
                    hold(p.viewObj.mainAxes, 'off')
                end
%                 maskIDs(maskIDs == 0) = [];
%                 for i = 1:numel(maskIDs)
%                     drawfreehand(p.viewObj.mainAxes, 'Position', maskTableTmp{maskTableTmp.maskID == maskIDs(i), {'y', 'x'}},...
%                         'Color', 'red', 'InteractionsAllowed', 'none');
%                 end
            else
                delete(p.imgMaskH)
                delete(p.cellMaskH)
            end
        end
        
        function scatterCallback(p, ~, ~)
            if p.viewObj.scatterCheckBox.Value == 0
                p.updateImageInView();
            end
            p.updateMainAxes();
        end
        
        function shuffleColorsInView(p, ~, ~)
            if ~logical(p.viewObj.scatterCheckBox.Value)
                nucTableTmp = p.IFboundaries.getNucBoundariesInRect(p.IFboundaries.channels(p.channelIdx), p.viewRect);
                cellTableTmp = p.IFboundaries.getCellBoundariesInRect(p.IFboundaries.channels(p.channelIdx), p.viewRect);
                tmpCellIDs = union(nucTableTmp.cellID, cellTableTmp.cellID);
                                
                colorPal  = single(d2utils.distinguishable_colors(50));
                newColors = colorPal(randi(50, numel(tmpCellIDs), 1),:);
                
                %Update nuclei boundaries wtih new colors
                nucIdx = ismember(p.IFboundaries.nucBoundaries2.cellID, tmpCellIDs);
                p.IFboundaries.nucBoundaries2.colors(nucIdx, :) = newColors;
                %Update cell boundaries wtih new colors
                cellIdx = ismember(p.IFboundaries.cellBoundaries2.cellID, tmpCellIDs);
                p.IFboundaries.cellBoundaries2.colors(cellIdx, :) = newColors;
                %Update plots
                delete(p.nucleiPlotH);
                delete(p.cellPlotH);
                p.overlayNuclei;
                p.overlayCells;
            end
        end
        
        function addNucAndCell(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucButton, p.viewObj.addCellButton,...
                p.viewObj.maskCellButton, p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.deleteNucButton,...
                p.viewObj.deleteCellButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            p.nucH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            if ~isempty(p.nucH.Position) %Allows 'escape' from ROI
                p.IFtable.addNucleus(fliplr(round(p.nucH.Position)), [0 0], true);
                p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
                p.updateCentroidListView();
                p.updateMainAxes();
            end
            delete(p.nucH)
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucButton, p.viewObj.addCellButton,...
                p.viewObj.maskCellButton, p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.deleteNucButton,...
                p.viewObj.deleteCellButton, p.viewObj.deleteMaskButton], 'Enable', 'on')            
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function addEmptyNuc(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addCellButton,...
                p.viewObj.maskCellButton, p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.deleteNucButton,...
                p.viewObj.deleteCellButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            p.nucH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            if ~isempty(p.nucH.Position) %Allows 'escape' from ROI
                p.IFtable.addNucleus(fliplr(round(p.nucH.Position)), [0 0], false);
                p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
                p.updateCentroidListView();
                p.updateMainAxes();
            end
            delete(p.nucH)
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addCellButton,...
                p.viewObj.maskCellButton, p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.deleteNucButton,...
                p.viewObj.deleteCellButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function addCell(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskCellButton, p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.deleteNucButton,...
                p.viewObj.deleteCellButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            p.cellH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            if ~isempty(p.cellH.Position) %Allows 'escape' from ROI
                p.IFtable.addCell(fliplr(round(p.cellH.Position)), [0 0]);
                p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
                p.updateCentroidListView();
                p.updateMainAxes();
            end
            delete(p.cellH)
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskCellButton, p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.deleteNucButton,...
                p.viewObj.deleteCellButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function deleteNuc(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskCellButton, p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.addCellButton,...
                p.viewObj.deleteCellButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            if p.deleteROI
                tmpMask = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
                if ~isempty(tmpMask.Position)
                    p.IFtable.deleteNuc(tmpMask.Position, p.viewRect, true);
                end
            else
                [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.
                if ~isempty(x)
                    p.IFtable.deleteNuc([x, y], p.viewRect, false);
                end
            end
            p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
            p.updateCentroidListView();
            p.updateMainAxes();
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskCellButton, p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.addCellButton,...
                p.viewObj.deleteCellButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function deleteCell(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskCellButton, p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.addCellButton,...
                p.viewObj.deleteNucButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            if p.deleteROI
                tmpMask = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
                if ~isempty(tmpMask.Position)
                    p.IFtable.deleteCell(tmpMask.Position, p.viewRect, true);
                end
            else
                [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.
                if ~isempty(x)
                    p.IFtable.deleteCell([x, y], p.viewRect, false);
                end
                p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
                p.updateCentroidListView();
                p.updateMainAxes();
            end
            p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
            p.updateCentroidListView();
            p.updateMainAxes();
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskCellButton, p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.addCellButton,...
                p.viewObj.deleteNucButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function deleteNucAndCell(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskCellButton, p.viewObj.maskImgButton, p.viewObj.deleteCellButton, p.viewObj.addCellButton,...
                p.viewObj.deleteNucButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            if p.deleteROI
                tmpMask = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
                if ~isempty(tmpMask.Position)
                    p.IFtable.deleteNucAndCell(tmpMask.Position, p.viewRect, true);
                end
            else
                [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.
                if ~isempty(x)
                    p.IFtable.deleteNucAndCell([x, y], p.viewRect, false);
                end
                p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
                p.updateCentroidListView();
                p.updateMainAxes();
            end
            p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
            p.updateCentroidListView();
            p.updateMainAxes();
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskCellButton, p.viewObj.maskImgButton, p.viewObj.deleteCellButton, p.viewObj.addCellButton,...
                p.viewObj.deleteNucButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function addImgMask(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskCellButton, p.viewObj.deleteNucAndCellButton, p.viewObj.deleteCellButton, p.viewObj.addCellButton,...
                p.viewObj.deleteNucButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            channel = p.IFtable.channels{p.channelIdx};
            p.maskH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            %addlistener(p.maskH, 'DrawingFinished', @p.maskSpots);
            if ~isempty(p.maskH.Position) && isvalid(p.maskH) %Allows 'escape' from ROI
                warning('off', 'MATLAB:polyshape:repairedBySimplify')
                tmpMaskPoly = polyshape(round(p.maskH.Position)); %This can reduce # of vertices. OK to delete. 
                p.maskImg(tmpMaskPoly.Vertices, channel)
                warning('on', 'MATLAB:polyshape:repairedBySimplify')
            end
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskCellButton, p.viewObj.deleteNucAndCellButton, p.viewObj.deleteCellButton, p.viewObj.addCellButton,...
                p.viewObj.deleteNucButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function maskImg(p, maskPosition, channel)
            %tmpPoly = roi.Position;
            set(p.viewObj.masksCheckBox, 'Value', true)
            p.IFtable.maskObj.addMaskLocalCoords(maskPosition, channel);
            delete(p.maskH)
            %Find cells and nuclei in maskBB
            tmpBB = d2utils.polygonBoundingBox(fliplr(maskPosition));
            p.IFtable.requantCellsInBB({tmpBB}, p.IFtable.channels(p.channelIdx));
            p.IFtable.updateCentroidList(p.IFtable.channels(p.channelIdx), p.quantMetricDict(p.quantMetric));
            p.updateCentroidListView();
            p.updateMainAxes();
        end
        
        function addCellMask(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.deleteCellButton, p.viewObj.addCellButton,...
                p.viewObj.deleteNucButton, p.viewObj.deleteMaskButton], 'Enable', 'off')
            channel = p.IFtable.channels{p.channelIdx};
            p.maskH = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
            if ~isempty(p.maskH.Position) && isvalid(p.maskH) %Allows 'escape' from ROI
                warning('off', 'MATLAB:polyshape:repairedBySimplify')
                tmpMaskPoly = polyshape(round(p.maskH.Position)); %This can reduce # of vertices. OK to delete. 
                p.maskCellInChannel(tmpMaskPoly.Vertices, channel)
                warning('on', 'MATLAB:polyshape:repairedBySimplify')
            end
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.deleteCellButton, p.viewObj.addCellButton,...
                p.viewObj.deleteNucButton, p.viewObj.deleteMaskButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function maskCellInChannel(p, maskPosition, channel)
            set(p.viewObj.masksCheckBox, 'Value', true)
            newMaskID = p.IFboundaries.maskObj.addMaskLocalCoords(maskPosition, channel);
            delete(p.maskH)
            cellIDs = p.IFboundaries.addCellMask(channel, maskPosition);
            p.IFtable.addCellMask(channel, cellIDs, newMaskID);
            p.IFtable.updateCentroidList(p.IFtable.channels(p.channelIdx), p.quantMetricDict(p.quantMetric));
            p.updateCentroidListView();
            p.updateMainAxes();
        end
        
        function deleteMask(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.deleteCellButton, p.viewObj.addCellButton,...
                p.viewObj.deleteNucButton, p.viewObj.maskCellButton], 'Enable', 'off')
            if p.deleteROI
                tmpMask = drawfreehand(p.viewObj.mainAxes, 'Parent', p.viewObj.mainAxes);
                tmpPoints = tmpMask.Position;
            else
                [x, y] = getpts(p.viewObj.mainAxes); %Simple but not interruptible. Can make WindowButtonDownFcn if want something interruptiblef.
                tmpPoints = d2utils.getPtsInsideView([x, y], p.viewRect);
            end
            if ~isempty(tmpPoints)
                [cellMaskIDs, ~, ~] = p.IFboundaries.maskObj.removeMasksByLocalPoints(tmpPoints, p.viewRect);
                if ~isempty(cellMaskIDs)
                    p.IFtable.removeCellMasks(cellMaskIDs, p.IFboundaries.channels(p.channelIdx));
                end
                [~, ~, imgMaskBB] = p.IFtable.maskObj.removeMasksByLocalPoints(tmpPoints, p.viewRect);
                if ~isempty(imgMaskBB)
                    p.IFtable.requantCellsInBB(imgMaskBB, p.IFtable.channels(p.channelIdx));
                end
                p.IFtable.updateCentroidList(p.IFtable.channels(p.channelIdx), p.quantMetricDict(p.quantMetric));
                p.updateCentroidListView();
                p.updateMainAxes();
            end
            set([p.viewObj.zoomAxes, p.viewObj.panAxes, p.viewObj.shuffleColors, p.viewObj.addNucAndCellButton, p.viewObj.addNucButton,...
                p.viewObj.maskImgButton, p.viewObj.deleteNucAndCellButton, p.viewObj.deleteCellButton, p.viewObj.addCellButton,...
                p.viewObj.deleteNucButton, p.viewObj.maskCellButton], 'Enable', 'on')
            set(p.viewObj.figHandle, 'WindowButtonDownFcn', {@p.figWindowDown})
        end
        
        function type = getSelectionType(p)
            type = get(p.viewObj.figHandle, 'SelectionType');
        end
        
        function p = changeCytoRadius(p, rad, varargin) %varargin for specifying channels
            %change radius
            p.IFtable.radius = rad;
            %requant cyto. 
            disp('Quantifying nuclei and cytoplasmic donut.')
            if nargin == 2
                p.IFtable.quantAllLabelMat2();
            elseif nargin ==3
                p.IFtable.quantAllLabelMat2(varargin{1});
            end
            %update centroid list
            cellIdx = get(p.viewObj.centroidList, 'Value');
            cellID = p.IFtable.centroidLists{p.channelIdx}.cellID(cellIdx);
            p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
            p.updateCentroidListView();
            newCellIdx = find(p.IFtable.centroidLists{p.channelIdx}.cellID == cellID);
            set(p.viewObj.centroidList, 'Value', newCellIdx);
            %update main axes
            p.updateMainAxes();
        end
        
        function p = adjustNucleiParameters(p, minSize, varargin)
            n = inputParser;
            n.addRequired('minSize', @isnumeric)
            n.addParameter('sensitivity', 0.1, @(x)validateattributes(x,{'numeric'}, {'scalar', '>=',0,'<=',1.0}));
            n.addParameter('strelRad', 40, @(x)validateattributes(x,{'numeric'}, {'scalar', 'nonnegative'}));
            n.addParameter('withCyto', true, @islogical);
            n.parse(minSize, varargin{:});
            %Set parameters
            p.IFboundaries.minNucArea = n.Results.minSize;
            p.IFboundaries.strelRadius = n.Results.strelRad;
            %Restitch dapi mask
            disp('masking dapi')
            if all(p.scanObj.scanDim == 1)
                p.IFboundaries.stitchDAPImask2('sensitivity', n.Results.sensitivity);
            else
                p.IFboundaries.stitchDAPImask(n.Results.sensitivity);
            end
            disp('Finding nuclei boundaries.')
            p.IFboundaries.makeNucleiLabelMat();
            
            if n.Results.withCyto
                if isfile(sprintf('%s_outlines.tif', p.IFboundaries.cellPoseCytoFile))
                    disp('Loading cellpose cytoplasmic boundaries.')
                    p.IFboundaries.loadCellPoseCyto(sprintf('%s_masks.tif', p.IFboundaries.cellPoseCytoFile), sprintf('%s_outlines.txt',  p.IFboundaries.cellPoseCytoFile));
                    p.IFboundaries.labelMat2cytoTable();
                    disp('Quantifying nuclei and cytoplasmic IF signal.')
                    p.IFtable.quantAllCytoBoundaries(true); %only quantifying nuclei overlapping cellpose cytoplasm masks
                else
                    disp('Quantifying nuclei and cytoplasmic donut.')
                    p.IFtable.quantAllLabelMat2();
                end
            else
                disp('Quantifying nuclei IF signal.')
                p.IFboundaries.addColors();
                p.IFtable.quantBoundaries();
            end
            %update centroid list
            cellIdx = get(p.viewObj.centroidList, 'Value');
            cellID = p.IFtable.centroidLists{p.channelIdx}.cellID(cellIdx);
            p.IFtable.makeCentroidList(p.quantMetricDict(p.quantMetric));
            p.updateCentroidListView();
            newCellIdx = find(p.IFtable.centroidLists{p.channelIdx}.cellID == cellID);
            set(p.viewObj.centroidList, 'Value', newCellIdx);
            %update main axes
            p.updateMainAxes();
        end
        
        function keyPressFunctions(p, ~, evt)
            keyPressed = evt.Key;
            modifierPressed = evt.Modifier;
            if isempty(modifierPressed)
                switch(keyPressed)
                    case 'z'
                        p.zoomInPressed();
                    case 'p'
                        p.panViewPressed();
                    case 'd'
                        checkBoxValue = get(p.viewObj.dapiCheckBox, 'Value');
                        set(p.viewObj.dapiCheckBox, 'Value', ~checkBoxValue);
                        p.updateMainAxes();
                    case 'f'
                        checkBoxValue = get(p.viewObj.IFCheckBox, 'Value');
                        set(p.viewObj.IFCheckBox, 'Value', ~checkBoxValue);
                        p.updateMainAxes();
                    case 'n'
                        checkBoxValue = get(p.viewObj.nucleiBordersCheckBox, 'Value');
                        set(p.viewObj.nucleiBordersCheckBox, 'Value', ~checkBoxValue);
                        p.overlayNuclei();
                    case 'c'
                        checkBoxValue = get(p.viewObj.cellBordersCheckBox, 'Value');
                        set(p.viewObj.cellBordersCheckBox, 'Value', ~checkBoxValue);
                        p.overlayCells();
                    case 's'
                        p.shuffleColorsInView();
                    case 'm'
                        p.addCellMask();
                    case 'i'
                        p.addImgMask();
                    case 'a'
                        p.addNucAndCell();
                    case 'w'
                        p.addEmptyNuc();
                    case 'e'
                        p.addCell();
                    case 'h'
                        p.deleteNucAndCell();
                    case 'j'
                        p.deleteNuc();
                    case 'k'
                        p.deleteCell();
                    case 'l'
                        p.deleteMask();
                    case 'g' 
                        p.fixedZoom = true;
                    case 'x'
                        set(p.viewObj.figHandle, 'WindowButtonDownFcn', '')
                    case 'uparrow'
                        cellIdx = max(1, get(p.viewObj.centroidList, 'Value')-1);
                        cellPos = p.IFtable.centroidLists{p.channelIdx}{cellIdx, {'x', 'y'}};
                        p.viewRect = d2utils.getRectAroundPoint(cellPos, 2 * p.cellViewRadius, 2 * p.cellViewRadius, p.scanObj.stitchDim);
                        set(p.viewObj.scatterCheckBox, 'Value', 0)
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.thumbCntrlr.overlayThumbnailRect();
                    case 'downarrow'
                        cellIdx = min(get(p.viewObj.centroidList, 'Value')+1, p.numCells);
                        cellPos = p.IFtable.centroidLists{p.channelIdx}{cellIdx, {'x', 'y'}};
                        p.viewRect = d2utils.getRectAroundPoint(cellPos, 2 * p.cellViewRadius, 2 * p.cellViewRadius, p.scanObj.stitchDim);
                        set(p.viewObj.scatterCheckBox, 'Value', 0)
                        p.updateImageInView();
                        p.updateMainAxes();
                        p.thumbCntrlr.overlayThumbnailRect();
                end
            elseif strcmp(modifierPressed{1}, 'shift')
                switch(keyPressed)
                    case 's'
                        p.viewObj.saveButtonPressed();
                end
            end
        end
    end
end