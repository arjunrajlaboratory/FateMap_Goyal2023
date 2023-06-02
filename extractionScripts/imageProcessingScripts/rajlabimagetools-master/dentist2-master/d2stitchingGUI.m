classdef d2stitchingGUI < handle
     
    properties (Access = public)
        %Model objects
        scanObj
        
        stitchCntrlr
        
        figHandle
        imageAxesTopL
        imageAxesBottomR
        imageAxesBottomL
        
        startPos
        startPosButtons
        snake
        snakeButtons
        direction
        directionButtons
        
        newPositionButton
        rowValue
        colValue
        rowText
        colText
        
        cpSelectRowButton
        cpSelectColButton
        previewStitchButton
    end
     
    methods
        
        function p = d2stitchingGUI(scanDim, scanFile) %Consider updating this with argument parser and default filename.
            validateattributes(scanDim, {'numeric'}, {'size', [1,2]})
            if isfile(scanFile)
                p.scanObj = scanObject('scanDim', scanDim, 'scanFile', scanFile);
            
                createComponents(p)
            
                startupFcn(p)
            else
                fprintf('Unable to find %s. Please input valid scan filename\n', scanFile)
                return
            end
        end
        
        function createComponents(p)
            p.figHandle = figure('Visible', 'off', 'Position', [100 100 800 800]);
            p.imageAxesBottomL = axes('Parent', p.figHandle, 'XTickLabel', '', 'YTickLabel', '', 'Position', [0.025 0.025 0.475 0.475], 'Interactions',[]);
            p.imageAxesTopL = axes('Parent', p.figHandle,'XTickLabel', '', 'YTickLabel', '', 'Position', [0.025 0.5 0.475 0.475], 'Interactions',[]);
            p.imageAxesBottomR = axes('Parent', p.figHandle,'XTickLabel', '', 'YTickLabel', '', 'Position', [0.5 0.025 0.475 0.475], 'Interactions',[]);
            
            p.startPos = uibuttongroup(p.figHandle, 'Units', 'normalized', 'Position', [0.525 0.85 0.15 0.1]);
            p.startPosButtons(1) = uicontrol(p.startPos, 'Style','radiobutton', 'String', 'top left', 'Units', 'normalized', 'Position', [0.01 0.8 0.9 0.2]);
            p.startPosButtons(2) = uicontrol(p.startPos, 'Style','radiobutton', 'String', 'top right' ,'Units', 'normalized', 'Position', [0.01 0.55 0.9 0.2]);
            p.startPosButtons(3) = uicontrol(p.startPos, 'Style','radiobutton', 'String', 'bottom left' ,'Units', 'normalized', 'Position', [0.01 0.3 0.9 0.2]);
            p.startPosButtons(4) = uicontrol(p.startPos, 'Style','radiobutton', 'String', 'bottom right' ,'Units', 'normalized', 'Position', [0.01 0.05 0.9 0.2]);
            
            p.snake = uibuttongroup(p.figHandle, 'Units', 'normalized', 'Position', [0.525 0.76 0.15 0.06]);
            p.snakeButtons(1) = uicontrol(p.snake, 'Style','radiobutton', 'String', 'snake', 'Units', 'normalized', 'Position', [0.01 0.5 0.9 0.4]);
            p.snakeButtons(2) = uicontrol(p.snake, 'Style','radiobutton', 'String', 'no snake' ,'Units', 'normalized', 'Position', [0.01 0.1 0.9 0.4]);
            
            p.direction = uibuttongroup(p.figHandle, 'Units', 'normalized', 'Position', [0.525 0.67 0.15 0.06]);
            p.directionButtons(1) = uicontrol(p.direction, 'Style','radiobutton', 'String', 'horizontal', 'Units', 'normalized', 'Position', [0.01 0.5 0.9 0.4]);
            p.directionButtons(2) = uicontrol(p.direction, 'Style','radiobutton', 'String', 'vertical' ,'Units', 'normalized', 'Position', [0.01 0.1 0.9 0.4]);
            
            p.newPositionButton = uicontrol('Style', 'pushbutton', 'String', 'show new positions', 'Units', 'normalized', 'Position', [0.69 0.90 0.14 0.05]);
            p.rowValue = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.735 0.86 0.08 0.035]);
            p.colValue = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.735 0.81 0.08 0.035]);
            p.rowText = uicontrol('Style', 'text', 'String', 'row', 'Units', 'normalized', 'Position', [0.69 0.85 0.04 0.035]);
            p.colText = uicontrol('Style', 'text', 'String', 'col', 'Units', 'normalized', 'Position', [0.69 0.8 0.04 0.035]);
            
            p.cpSelectRowButton = uicontrol('Style', 'pushbutton', 'String', ['<html><center>select row<br />control points</center><html>'], 'Units', 'normalized', 'Position', [0.84 0.90 0.14 0.05]);
            p.cpSelectColButton = uicontrol('Style', 'pushbutton', 'String', ['<html><center>select col<br />control points</center><html>'], 'Units', 'normalized', 'Position', [0.84 0.85 0.14 0.05]);
            p.previewStitchButton = uicontrol('Style', 'pushbutton', 'String', 'preview stitch', 'Units', 'normalized', 'Position', [0.84 0.8 0.14 0.05]);

            p.figHandle.Visible = 'on';
            
        end
        
        function startupFcn(p)
            p.stitchCntrlr = d2stitchController(p, p.scanObj);  
            p.attachController(p.stitchCntrlr);
        end
        
        function attachController(p, controller)
            p.startPos.SelectionChangedFcn = {@controller.startSelectionChanged};
            p.snake.SelectionChangedFcn = {@controller.snakeSelectionChanged};
            p.direction.SelectionChangedFcn = {@controller.directionSelectionChanged};
            p.newPositionButton.Callback = {@controller.newPositionPushed};
            p.rowValue.Callback = {@controller.updateImagesInView};
            p.colValue.Callback = {@controller.updateImagesInView};
            p.cpSelectRowButton.Callback = {@controller.cpSelectRowPushed};
            p.cpSelectColButton.Callback = {@controller.cpSelectColPushed};
            p.previewStitchButton.Callback = {@controller.previewStitchPushed};
            p.figHandle.CloseRequestFcn = {@p.closeFigFcn};
        end
        
        function closeFigFcn(p, ~, ~)
            delete(p.figHandle)
            if isempty(p.stitchCntrlr.rowTransform)
                p.scanObj.rowTransformCoords = [0, round(0.9 * p.scanObj.tileSize(1))]; %Set default overlap to 10%
            else
                p.scanObj.rowTransformCoords = round(median(p.stitchCntrlr.rowTransform, 1));
            end
            
            if isempty(p.stitchCntrlr.colTransform)
                p.scanObj.columnTransformCoords = [round(0.9 * p.scanObj.tileSize(2)), 0]; %Set default overlap to 10%
            else
                p.scanObj.columnTransformCoords = round(median(p.stitchCntrlr.colTransform, 1));
            end
            p.scanObj.snake = p.stitchCntrlr.snakeValue;
            p.scanObj.startPos = p.stitchCntrlr.startPosValue;
            p.scanObj.direction = p.stitchCntrlr.directionValue;
            p.scanObj.loadTiles();
            disp('saving tilesTable.csv')
            p.scanObj.saveTilesTable();
            
            if isfile('scanSummary.txt') %In case we want to update scanSummary with new transformCoords but don't want to change spot thresholds
                disp('updating scanSummary.txt')
                scanSummaryTable = d2utils.parseScanSummary('scanSummary.txt');
                scanSummaryTable{'stitchDimensions',:} = {num2str(p.scanObj.stitchDim)};
                scanSummaryTable{'startPosition',:} = {p.scanObj.startPos};
                scanSummaryTable{'scanDirection',:} = {p.scanObj.direction};
                scanSummaryTable{'snake',:} = {string(p.scanObj.snake)};
                scanSummaryTable{'rowTransform',:} = {num2str(p.scanObj.rowTransformCoords)};
                scanSummaryTable{'columnTransform',:} = {num2str(p.scanObj.columnTransformCoords)};
                writetable(scanSummaryTable, 'scanSummary.txt', 'WriteRowNames', true, 'WriteVariableNames', false, 'QuoteStrings', false, 'Delimiter', '\t') 
            else
                disp('saving scanSummary.txt')
                p.scanObj.saveScanSummary();
            end
        end
    
    end
end