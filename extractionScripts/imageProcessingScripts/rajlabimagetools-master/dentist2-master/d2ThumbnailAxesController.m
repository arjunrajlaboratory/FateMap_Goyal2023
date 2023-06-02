classdef d2ThumbnailAxesController < handle
    
    properties
        
        %Model objects
        spotTable
        scanObj
        
        %View
        viewObj
        
        %Other controllers
        mainAxesCntrlr
        threshCntrlr
        
        thumbPlotH
        thumbRectH
        
    end
    
    methods
        function p = d2ThumbnailAxesController(view, scanObj, spotTable)
            p.viewObj = view;
            p.spotTable = spotTable;
            p.scanObj = scanObj;
        end
        
        function startup(p)
            p.plotThumbnail();
        end
        
        function plotThumbnail(p)
            set(p.viewObj.thumbAxes, 'XLim', [1,  p.scanObj.stitchDim(2)])
            set(p.viewObj.thumbAxes, 'YLim', [1,  p.scanObj.stitchDim(1)])
            hold(p.viewObj.thumbAxes, 'on')
            %Note binning about every 1000 pixels. Can change this with the
            %third argument to binscatter
            p.thumbPlotH = binscatter(p.spotTable.centroidLists{p.mainAxesCntrlr.channelIdx}.y, p.spotTable.centroidLists{p.mainAxesCntrlr.channelIdx}.x, round(min(p.scanObj.stitchDim/1000)),...
                'Parent', p.viewObj.thumbAxes, 'HitTest','off');
            p.overlayThumbnailRect();
            hold(p.viewObj.thumbAxes, 'off')
            
        end
        
        function overlayThumbnailRect(p)
            if any(p.mainAxesCntrlr.viewRect(3:4) < p.scanObj.stitchDim) && all(p.mainAxesCntrlr.viewRect(3:4) > 4) %No overlay if too zoomed out or in
                if isempty(p.thumbRectH) || ~isvalid(p.thumbRectH)
                    p.thumbRectH = rectangle(...
                    'Position', d2utils.rotateRectROI(p.mainAxesCntrlr.viewRect), ...
                    'EdgeColor', 'b', 'LineWidth', 2,...
                    'Parent', p.viewObj.thumbAxes,...
                    'Hittest', 'off');
                else
                    set(p.thumbRectH, 'Position', d2utils.rotateRectROI(p.mainAxesCntrlr.viewRect))
                end
            else
                delete(p.thumbRectH)
            end
        end
        
        function thumbAxesButtonDown(p, ~, ~)
            if any(p.mainAxesCntrlr.viewRect(3:4) < p.scanObj.stitchDim)
                currentPoint = get(p.viewObj.thumbAxes, 'CurrentPoint');
                switch(p.mainAxesCntrlr.getSelectionType)
                    case 'open'
                        newRect = d2utils.getRectAroundPoint(currentPoint(1,2:-1:1), p.mainAxesCntrlr.viewRect(3), p.mainAxesCntrlr.viewRect(4), p.scanObj.stitchDim);
                        p.mainAxesCntrlr.viewRect = newRect;
                        p.mainAxesCntrlr.updateImageInView();
                        p.mainAxesCntrlr.updateMainAxes();
                        p.overlayThumbnailRect();
                    case 'normal'
                        if d2utils.pointInSideViewRect(p.mainAxesCntrlr.viewRect, currentPoint)
                            set(p.viewObj.figHandle, 'WindowButtonUpFcn', {@p.stopThumbDrag});
                            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', {@p.dragThumb});
                        end
                end
            end
        end
        
        function dragThumb(p, ~, ~)
            currentPoint = get(p.viewObj.thumbAxes, 'CurrentPoint');
            newRect = d2utils.getRectAroundPoint(currentPoint(1,2:-1:1), p.mainAxesCntrlr.viewRect(3), p.mainAxesCntrlr.viewRect(4), p.scanObj.stitchDim);
            p.mainAxesCntrlr.viewRect = newRect;
            p.overlayThumbnailRect();
            p.mainAxesCntrlr.updateImageInView();
            p.mainAxesCntrlr.updateMainAxes();
        end
        
        function stopThumbDrag(p, ~, ~)
            set(p.viewObj.figHandle, 'WindowButtonMotionFcn', '');
            set(p.viewObj.figHandle, 'WindowButtonUpFcn', '');
        end
    end
end