classdef d2stitchController < handle
    
    properties
        
        scanObj
        viewObj
        
        startPosValue
        snakeValue
        directionValue
        
        rowTransform = []
        colTransform = []
    end
    methods
        function p = d2stitchController(view, scanObject)
            p.viewObj = view;
            p.scanObj = scanObject;
            if ~isempty(p.scanObj.rowTransformCoords)
                p.rowTransform = p.scanObj.rowTransformCoords;
            end
            
            if ~isempty(p.scanObj.columnTransformCoords)
                p.colTransform = p.scanObj.columnTransformCoords;
            end
            p.startup()
   
        end
        
        function startup(p)
            set(p.viewObj.rowValue, 'String', num2str(max(2,floor(p.scanObj.scanDim(1)/2)))) %Set default view to middle of scan, with a minimum of 2 to avoid an index-1 issue later on
            set(p.viewObj.colValue, 'String', num2str(max(2,floor(p.scanObj.scanDim(2)/2))))
            p.startPosValue = get(p.viewObj.startPos.SelectedObject, 'String');
            p.snakeValue = strcmp(get(p.viewObj.snake.SelectedObject, 'String'), 'snake');
            p.directionValue = get(p.viewObj.direction.SelectedObject, 'String');
            p.scanObj.scanMatrix = d2utils.makeScanMatrix(p.scanObj.scanDim, 'start', p.startPosValue, 'snake', p.snakeValue,'direction', p.directionValue); 
%             p.scanObj.scanMatrix = d2utils.makeScanMatrix(p.scanObj.scanDim);            
            p.updateImagesInView();
        end
        
        function updateImagesInView(p, ~, ~)
            row =  str2double(get(p.viewObj.rowValue, 'String'));
            col = str2double(get(p.viewObj.colValue, 'String'));
            tileLL = p.scanObj.scanMatrix(row, col);
            tileUL = p.scanObj.scanMatrix(row-1, col);
            tileLR = p.scanObj.scanMatrix(row, col+1);
            
            imgLL = p.scanObj.getTileFromScan(tileLL, 'dapi');
            imgUL = p.scanObj.getTileFromScan(tileUL, 'dapi');
            imgLR =  p.scanObj.getTileFromScan(tileLR, 'dapi');
            %If desired, could add contrasting control here.  
            imshow(scale(imgLL), 'Parent', p.viewObj.imageAxesBottomL)
            imshow(scale(imgUL), 'Parent', p.viewObj.imageAxesTopL)
            imshow(scale(imgLR), 'Parent', p.viewObj.imageAxesBottomR)
        end
        
        function startSelectionChanged(p, ~, evt)
            p.startPosValue = evt.NewValue.String;
            p.scanObj.scanMatrix = d2utils.makeScanMatrix(p.scanObj.scanDim, 'start', p.startPosValue, 'snake', p.snakeValue,'direction', p.directionValue); 
            p.updateImagesInView();
        end
        
        function snakeSelectionChanged(p, ~, evt)
            p.snakeValue = strcmp(evt.NewValue.String, 'snake');
            p.scanObj.scanMatrix = d2utils.makeScanMatrix(p.scanObj.scanDim, 'start', p.startPosValue, 'snake', p.snakeValue,'direction', p.directionValue); 
            p.updateImagesInView();
        end
        
        function directionSelectionChanged(p, ~, evt)
            p.directionValue = evt.NewValue.String;
            p.scanObj.scanMatrix = d2utils.makeScanMatrix(p.scanObj.scanDim, 'start', p.startPosValue, 'snake', p.snakeValue,'direction', p.directionValue); 
            p.updateImagesInView();
        end
        
        function newPositionPushed(p, ~, ~)
            set(p.viewObj.rowValue, 'String', int2str(randi([2, p.scanObj.scanDim(1)])))
            set(p.viewObj.colValue, 'String', int2str(randi(p.scanObj.scanDim(2)-1)))
            p.updateImagesInView();
        end
        
        function cpSelectRowPushed(p, ~, ~)
            row =  str2double(get(p.viewObj.rowValue, 'String'));
            col = str2double(get(p.viewObj.colValue, 'String'));
            tileLL = p.scanObj.scanMatrix(row, col);
            tileUL = p.scanObj.scanMatrix(row-1, col);
            imgLL = p.scanObj.getTileFromScan(tileLL, 'dapi');
            imgUL = p.scanObj.getTileFromScan(tileUL, 'dapi');
            [moving_out,fixed_out] = cpselect(scale(imgLL),scale(imgUL),'Wait',true);
            if ~isempty(moving_out)
                %Not sure if it's worth using cpcorr here to try improving
                %control points. 
                p.rowTransform = [p.rowTransform; median(fixed_out-moving_out, 1)];
            end
        end
        
        function cpSelectColPushed(p, ~, ~)
            row =  str2double(get(p.viewObj.rowValue, 'String'));
            col = str2double(get(p.viewObj.colValue, 'String'));
            tileLL = p.scanObj.scanMatrix(row, col);
            tileLR = p.scanObj.scanMatrix(row, col+1);
            imgLL = p.scanObj.getTileFromScan(tileLL, 'dapi');
            imgLR =  p.scanObj.getTileFromScan(tileLR, 'dapi');
            [moving_out,fixed_out] = cpselect(scale(imgLL),scale(imgLR),'Wait',true);
            if ~isempty(moving_out)
                %Not sure if it's worth using cpcorr here to try improving
                %control points. 
                p.colTransform = [p.colTransform; median(moving_out-fixed_out, 1)];
            end
        end
        
        function previewStitchPushed(p, ~, ~)
            row =  str2double(get(p.viewObj.rowValue, 'String'));
            col = str2double(get(p.viewObj.colValue, 'String'));
            if isempty(p.rowTransform)
                p.rowTransform = [0, 0.9 * p.scanObj.tileSize(1)]; %Set default overlap to 10%
            end
            
            if isempty(p.colTransform)
                p.colTransform = [0.9 * p.scanObj.tileSize(2), 0]; %Set default overlap to 10%
            end
            
            tmpStitch = p.scanObj.stitchTiles([row-1,row], [col,col+1], 'dapi', round(median(p.rowTransform, 1)), round(median(p.colTransform, 1)));
            disp(size(tmpStitch))
            figure
            imshow(scale(tmpStitch));
        end
    end
    
end
