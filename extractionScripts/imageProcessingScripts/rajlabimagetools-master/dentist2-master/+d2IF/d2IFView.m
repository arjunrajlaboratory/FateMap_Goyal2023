classdef d2IFView < handle
    
    properties (Access = public)
        
        %Model objects
        scanObj
        maskObj
        IFboundaries
        IFtable
        
        mainAxesCntrlr
        thumbAxesCntrlr
        
        figHandle
        mainAxes 
        thumbAxes
        centroidList %List box 
        quantMetric
        quantLabel
        quantButtons
        selectionTool
        selectionButton
        selectionLabel
        nucleiBordersCheckBox %Check box
        cellBordersCheckBox %Check box
        scatterCheckBox %Check box
        masksCheckBox  %Check box
        dapiCheckBox  %Check box
        IFCheckBox  %Check box
        
        channelPopup %Drop down
        colormapPopup %Drop down
        saveButton %Push button
        exportButton %Push button
        zoomAxes %Push button
        panAxes %Push button
        shuffleColors %Push button
        
        addNucButton %Push button
        addNucAndCellButton %Push button
        addCellButton %Push button
        deleteNucButton %Push button
        deleteCellButton %Push button
        deleteNucAndCellButton %Push button
        maskCellButton %Push button
        maskImgButton %Push button
        deleteMaskButton %Push button
        
        upperContrastSlider %Slider
        lowerContrastSlider %Slider
        sliderLabel
        
        zoomH
        zoomRect
        
    end
    
    % GUI startup and deletion
    methods (Access = public)

        % Construct app
        function p = d2IFView(scanObject, IFboundaries, IFtable)
            p.scanObj = scanObject;
            p.IFboundaries = IFboundaries;
            p.IFtable = IFtable;

            % Create UIFigure and components
            createComponents(p)
            
            startupFcn(p)

            if nargout == 0
                clear p
            end
        end

%         % Code that executes before app deletion
%         function delete(p)
%             %save tables
% 
%             % Delete UIFigure when app is deleted
%             delete(app.UIFigure)
%         end
    end
    
    % GUI construction
    methods 
        
        function createComponents(p)

            p.figHandle = figure('Visible', 'off', 'Position', [100 100 1550 900]);
            p.mainAxes = axes('Parent', p.figHandle, 'Position', [0.025 0.025 0.60 0.95], 'Ydir','reverse', 'XAxisLocation', 'bottom', 'YAxisLocation', 'left', 'Interactions',[]);
            set(p.mainAxes.Toolbar, 'Visible','off');
            p.thumbAxes = axes('Parent', p.figHandle, 'XTickLabel', '', 'YTickLabel', '', 'Ydir', 'reverse', 'Position', [0.64 0.645 0.21 0.33], 'Interactions', []);
            set(p.thumbAxes.Toolbar, 'Visible', 'off')
            set(p.thumbAxes, 'Colormap', flipud(hot))
            
            p.centroidList = uicontrol('Style', 'listbox', 'String', string(p.IFtable.centroidLists{1}{:,4}), 'Units', 'normalized', 'Position', [0.86 0.645 0.12 0.33]);
            p.quantMetric = uibuttongroup(p.figHandle, 'Units', 'normalized', 'Position', [0.86 0.51 0.12 0.1]);
            p.quantButtons(1) = uicontrol(p.quantMetric, 'Style','radiobutton', 'String', 'mean nucleus', 'Units', 'normalized', 'Position', [0.01 0.8 0.9 0.2]);
            p.quantButtons(2) = uicontrol(p.quantMetric, 'Style','radiobutton', 'String', 'mean cytoplasm' ,'Units', 'normalized', 'Position', [0.01 0.55 0.9 0.2]);
            p.quantButtons(3) = uicontrol(p.quantMetric, 'Style','radiobutton', 'String', 'sum nucleus' ,'Units', 'normalized', 'Position', [0.01 0.3 0.9 0.2]);
            p.quantButtons(4) = uicontrol(p.quantMetric, 'Style','radiobutton', 'String', 'sum cytoplasm' ,'Units', 'normalized', 'Position', [0.01 0.05 0.9 0.2]);
            p.quantLabel = uicontrol('Style', 'text', 'String', 'Quantification metric:', 'Units', 'normalized', 'Position', [0.86 0.61 0.07 0.02]);
            
            p.selectionTool = uibuttongroup(p.figHandle, 'Units', 'normalized', 'BorderType', 'none', 'Position', [0.87 0.45 0.1111 0.03]);
            p.selectionButton(1) = uicontrol(p.selectionTool, 'Style','radiobutton', 'String', 'points', 'Units', 'normalized', 'Position', [0.05 0.01 0.4 0.95]);
            p.selectionButton(2) = uicontrol(p.selectionTool, 'Style','radiobutton', 'String', 'draw' ,'Units', 'normalized', 'Position', [0.55 0.01 0.4 0.95]);
            p.selectionLabel = uicontrol('Style', 'text', 'String', 'Selection tool:', 'Units', 'normalized', 'Position', [0.885 0.475 0.07 0.02]);

            p.nucleiBordersCheckBox = uicontrol('Style', 'checkbox', 'String', 'nuclei (n)', 'Value', true, 'Units', 'normalized', 'Position', [0.65 0.604 0.0722 0.0333]);
            p.cellBordersCheckBox = uicontrol('Style', 'checkbox', 'String', 'cells (c)', 'Value', true, 'Units', 'normalized', 'Position', [0.71 0.604 0.0722 0.0333]);
            p.masksCheckBox = uicontrol('Style', 'checkbox', 'String', 'masks', 'Value', true, 'Units', 'normalized', 'Position', [0.76 0.604 0.0722 0.0333]);
            p.dapiCheckBox = uicontrol('Style', 'checkbox', 'String', 'DAPI (d)', 'Value', true, 'Units', 'normalized', 'Position', [0.65 0.58 0.0722 0.0333]);
            p.IFCheckBox = uicontrol('Style', 'checkbox', 'String', sprintf('%s (f)', p.IFtable.channels{1}),  'Value', true, 'Units', 'normalized', 'Position', [0.71 0.58 0.0722 0.0333]);
            p.scatterCheckBox = uicontrol('Style', 'checkbox', 'String', 'scatter', 'Value', true, 'Units', 'normalized', 'Position', [0.76 0.58 0.0722 0.0333]);

            p.upperContrastSlider = uicontrol('Style', 'slider', 'Value', 0.5, 'Units', 'normalized', 'Position', [0.65 0.56 0.1233 0.0050]);
            p.lowerContrastSlider = uicontrol('Style', 'slider', 'Units', 'normalized', 'Position', [0.65 0.53 0.1233 0.0050]);
            p.sliderLabel = uicontrol('Style', 'text', 'String', 'Adjust contrast', 'Units', 'normalized', 'Position', [0.680 0.48 0.07 0.03]);

            p.channelPopup = uicontrol('Style', 'popupmenu', 'String', p.IFtable.channels, 'Units', 'normalized', 'Position', [0.64 0.44 0.1111 0.0367]);
            p.colormapPopup = uicontrol('Style', 'popupmenu', 'String', p.IFtable.expressionColorPal, 'Units', 'normalized', 'Position', [0.755 0.44 0.1111 0.0367]);
            
            p.zoomAxes = uicontrol('Style', 'pushbutton', 'String', 'zoom (z)', 'Units', 'normalized', 'Position', [0.64 0.41 0.1111 0.0367]);
            p.panAxes = uicontrol('Style', 'pushbutton', 'String', 'pan view (p)', 'Units', 'normalized', 'Position', [0.64 0.37 0.1111 0.0367]);
            p.shuffleColors = uicontrol('Style', 'pushbutton', 'String', 'shuffle colors (s)', 'Units', 'normalized', 'Position', [0.64 0.33 0.1111 0.0367]);
            p.saveButton = uicontrol('Style', 'pushbutton', 'String', 'save/export (S)', 'Units', 'normalized', 'Position', [0.64 0.29 0.1111 0.0367]);

            p.addNucAndCellButton = uicontrol('Style', 'pushbutton', 'String', 'add nucleus & cell (a)', 'Units', 'normalized', 'Position', [0.755 0.41 0.1111 0.0367]);
            p.addNucButton = uicontrol('Style', 'pushbutton', 'String', 'add nucleus (w)', 'Units', 'normalized', 'Position', [0.755 0.37 0.1111 0.0367]);
            p.addCellButton = uicontrol('Style', 'pushbutton', 'String', 'add cell (e)', 'Units', 'normalized', 'Position', [0.755 0.33 0.1111 0.0367]);
            p.maskCellButton = uicontrol('Style', 'pushbutton', 'String', 'mask cells (m)', 'Units', 'normalized', 'Position', [0.755 0.29 0.1111 0.0367]);
            p.maskImgButton = uicontrol('Style', 'pushbutton', 'String', 'mask image (i)', 'Units', 'normalized', 'Position', [0.755 0.25 0.1111 0.0367]);

            p.deleteNucAndCellButton = uicontrol('Style', 'pushbutton', 'String', 'delete nuclei & cells (h)', 'Units', 'normalized', 'Position', [0.87 0.41 0.1111 0.0367]);
            p.deleteNucButton = uicontrol('Style', 'pushbutton', 'String', 'delete nuclei (j)', 'Units', 'normalized', 'Position', [0.87 0.37 0.1111 0.0367]);
            p.deleteCellButton = uicontrol('Style', 'pushbutton', 'String', 'delete cells (k)', 'Units', 'normalized', 'Position', [0.87 0.33 0.1111 0.0367]);
            p.deleteMaskButton = uicontrol('Style', 'pushbutton', 'String', 'delete masks (l)', 'Units', 'normalized', 'Position', [0.87 0.29 0.1111 0.0367]); 
%             p.exportButton = uicontrol('Style', 'pushbutton', 'String', 'export (E)', 'Units', 'normalized', 'Position', [0.87 0.41 0.1111 0.0367]);
            p.figHandle.Visible = 'on';
        end
        
        %attach controller
        function startupFcn(p)

            p.mainAxesCntrlr = d2IF.d2IFController(p, p.scanObj, p.IFboundaries, p.IFtable);  
            p.attatchMainAxesController(p.mainAxesCntrlr);        
            
            p.thumbAxesCntrlr = d2ThumbnailAxesController(p, p.scanObj, p.IFtable);  
            p.attachThumbnailController(p.thumbAxesCntrlr);
            
             %Connect controllers
            p.mainAxesCntrlr.thumbCntrlr = p.thumbAxesCntrlr;
            
            p.thumbAxesCntrlr.mainAxesCntrlr = p.mainAxesCntrlr;
            p.thumbAxesCntrlr.startup();
        end
        
        function p = attatchMainAxesController(p, controller)
            p.channelPopup.Callback = {@controller.changeChannel};
            p.colormapPopup.Callback = {@controller.changeColormap};
            p.centroidList.Callback = {@controller.centroidSelected};
            p.nucleiBordersCheckBox.Callback = {@controller.overlayNuclei};
            p.cellBordersCheckBox.Callback = {@controller.overlayCells};
            p.scatterCheckBox.Callback = {@controller.scatterCallback};
            p.dapiCheckBox.Callback = {@controller.updateMainAxes};
            p.IFCheckBox.Callback = {@controller.updateMainAxes};
            p.quantMetric.SelectionChangedFcn = {@controller.quantMetricChanged};
            p.selectionTool.SelectionChangedFcn = {@controller.selectionToolChanged};
            p.masksCheckBox.Callback = {@controller.overlayMasks};
            p.addNucAndCellButton.Callback = {@controller.addNucAndCell};
            p.addNucButton.Callback = {@controller.addEmptyNuc};
            p.addCellButton.Callback = {@controller.addCell};
            p.deleteNucAndCellButton.Callback = {@controller.deleteNucAndCell};
            p.deleteCellButton.Callback = {@controller.deleteCell};
            p.deleteNucButton.Callback = {@controller.deleteNuc};
            p.deleteMaskButton.Callback = {@controller.deleteMask};
            p.maskCellButton.Callback = {@controller.addCellMask};
            p.maskImgButton.Callback = {@controller.addImgMask};
            p.zoomAxes.Callback = {@controller.zoomInPressed};
            p.panAxes.Callback = {@controller.panViewPressed};
            p.upperContrastSlider.Callback = {@controller.updateMainAxes};
            p.lowerContrastSlider.Callback = {@controller.updateMainAxes};
            p.saveButton.Callback = {@p.saveButtonPressed};
            p.shuffleColors.Callback = {@controller.shuffleColorsInView};
%             
            p.figHandle.WindowButtonDownFcn = {@controller.figWindowDown};
            p.figHandle.KeyPressFcn = {@controller.keyPressFunctions};
%             p.figHandle.CloseRequestFcn = {@p.closeFigFcn}; 
            
            %Set KeyPressFcn for all uicontrols. 
%             set(findobj(p.figHandle, 'Type', 'UIControl'), 'KeyPressFcn', {@controller.keyPressFunctions})
        end
        
        function p = attachThumbnailController(p, controller)
            p.thumbAxes.ButtonDownFcn = {@controller.thumbAxesButtonDown};
        end
        
        function saveButtonPressed(p, ~, ~)
            fprintf('Saving IF table, IF boundaries, and masks.\nThis may take a minute...')
            p.IFtable.saveTable;
            p.IFtable.maskObj.saveMasksTable;
            p.IFboundaries.saveBoundaries;
            p.IFboundaries.maskObj.saveMasksTable;
            fprintf('done.\n')
        end

        function relaunchGUI(p)
            createComponents(p)
            startupFcn(p)
        end
        
%         function closeFigFcn(p, ~, ~)
%             delete(p.figHandle)
%             fprintf('Saving mask table, cell table, and spot tables.\nThis may take a minute\n')
%             p.saveButtonPressed;
%         end
        
    end

end
