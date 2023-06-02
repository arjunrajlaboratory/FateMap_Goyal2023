classdef d2ThresholdView2 < handle
    
    properties (Access = public)
        
        %Model objects
        spotTable
        scanObj
        maskObj
        nucleiObj
        
        mainAxesCntrlr
        threshAxesCntrlr
        thumbAxesCntrlr
        
        figHandle
        mainAxes 
        thumbAxes
        threshAxes
        centroidList %List box 
        spotsCheckBox %Check box
        centroidsCheckBox %Check box
        scatterCheckBox %Check box
        masksCheckBox  %Check box
        
        channelPopup %Drop down
        colormapPopup %Drop down
        saveButton %Push button
        exportButton %Push button
        zoomAxes %Push button
        panAxes %Push button
        shuffleColors %Push button
        
        zoomThresh %Push button
        filterMasksThresh %Push button
        threshValue %Edit field
        
        addCellButton %Toggle button
        deleteCellButton %Toggle button
        maskCellButton %Push button
        maskSpotsButton %Push button
        maskSpotsAllChannelsCheckBox % Check box
        maskSpotsForAllChannels = false
        deleteMaskButton %Push button
        
        upperContrastSlider %Slider
        lowerContrastSlider %Slider
        sliderLabel
        
        zoomH
        zoomRect
        
        showSpots = true %Not sure if these handles are necessary. 
        showCentroids = true
        showScatter= true
        showMasks = true
        
    end
    
    % GUI startup and deletion
    methods (Access = public)

        % Construct app
        function p = d2ThresholdView2(scanObject, maskObject, nucleiObject, spotTable)
            p.scanObj = scanObject;
            p.maskObj = maskObject;
            p.nucleiObj = nucleiObject;
            p.spotTable = spotTable;

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
            p.threshAxes = axes('Parent', p.figHandle, 'Position', [0.64 0.025 0.34 0.34], 'Interactions', []);
            set(p.threshAxes.Toolbar, 'Visible', 'off')
            
            %p.channelPopup = uicontrol('Style', 'popupmenu', 'String', p.spotTable.spotChannels, 'Units', 'normalized', 'Position', [0.64 0.49 0.1111 0.0367]); % CHANGE
            p.channelPopup = uicontrol('Style', 'popupmenu', 'String', p.scanObj.channels(~ismember(p.scanObj.channelTypes,'dapi')), 'Units', 'normalized', 'Position', [0.64 0.51 0.1111 0.0367]); % CHANGE
            p.colormapPopup = uicontrol('Style', 'popupmenu', 'String', p.spotTable.expressionColorPal, 'Units', 'normalized', 'Position', [0.64 0.48 0.1111 0.0367]);
            p.centroidList = uicontrol('Style', 'listbox', 'String', string(p.spotTable.centroidLists{1}.GroupCount),'Units', 'normalized', 'Position', [0.86 0.645 0.12 0.33]);
            p.spotsCheckBox = uicontrol('Style', 'checkbox', 'String', 'spots (s)', 'Value', p.showSpots,'Units', 'normalized', 'Position', [0.65 0.585 0.0722 0.0333]);
            p.centroidsCheckBox = uicontrol('Style', 'checkbox', 'String', 'nuclei (n)', 'Value', p.showCentroids, 'Units', 'normalized', 'Position', [0.695 0.585 0.0722 0.0333]);
            p.masksCheckBox = uicontrol('Style', 'checkbox', 'String', 'masks', 'Value', p.showMasks, 'Units', 'normalized', 'Position', [0.74 0.585 0.0722 0.0333]);
            p.scatterCheckBox = uicontrol('Style', 'checkbox', 'String', 'scatter (c)', 'Value', p.showScatter, 'Units', 'normalized', 'Position', [0.785 0.585 0.0722 0.0333]);
            
            %In case there are non-FISH channels, set channelPopup value
            %and centroidList to correspond to first FISH channel
            firstFISHidx = find(strcmp(p.channelPopup.String, p.spotTable.spotChannels{1}));
            set(p.channelPopup, 'Value', firstFISHidx)

            p.addCellButton = uicontrol('Style', 'pushbutton', 'String', 'add cells', 'Units', 'normalized', 'Position', [0.87 0.49 0.1111 0.0367]);
            p.deleteCellButton = uicontrol('Style', 'pushbutton', 'String', 'delete cells', 'Units', 'normalized', 'Position', [0.87 0.45 0.1111 0.0367]);

            p.maskSpotsAllChannelsCheckBox = uicontrol('Style', 'checkbox', 'String', sprintf('Mask all channels (a)'), 'Value', p.maskSpotsForAllChannels, 'Units', 'normalized', 'Position', [0.755 0.49 0.1111 0.0367]);
            p.maskSpotsButton = uicontrol('Style', 'pushbutton', 'String', 'mask spots (w)', 'Units', 'normalized', 'Position', [0.755 0.45 0.1111 0.0367]);
%             p.maskSpotsAllChannelsCheckBox = uicontrol('Style', 'pushbutton', 'String', 'mask spots all channels (a)', 'Units', 'normalized', 'Position', [0.755 0.45 0.1111 0.0367]);
            
            p.maskCellButton = uicontrol('Style', 'pushbutton', 'String', 'mask cells (x)', 'Units', 'normalized', 'Position', [0.755 0.41 0.1111 0.0367]);
            p.deleteMaskButton = uicontrol('Style', 'pushbutton', 'String', 'delete mask (d)', 'Units', 'normalized', 'Position', [0.755 0.37 0.1111 0.0367]); 

            p.zoomAxes = uicontrol('Style', 'pushbutton', 'String', 'zoom (z)', 'Units', 'normalized', 'Position', [0.64 0.45 0.1111 0.0367]);
            p.panAxes = uicontrol('Style', 'pushbutton', 'String', 'pan view (p)', 'Units', 'normalized', 'Position', [0.64 0.41 0.1111 0.0367]);
            p.saveButton = uicontrol('Style', 'pushbutton', 'String', 'save (S)', 'Units', 'normalized', 'Position', [0.87 0.41 0.1111 0.0367]);
            p.exportButton = uicontrol('Style', 'pushbutton', 'String', 'export (E)', 'Units', 'normalized', 'Position', [0.87 0.37 0.1111 0.0367]);
%             p.shuffleColors = uicontrol('Style', 'pushbutton', 'String', 'shuffle colors', 'Units', 'normalized', 'Position', [0.87 0.49 0.1111 0.0367]);
            p.shuffleColors = uicontrol('Style', 'pushbutton', 'String', 'shuffle colors', 'Units', 'normalized', 'Position', [0.64 0.37 0.1111  0.0367]);
            
            p.zoomThresh = uicontrol('Style', 'pushbutton', 'String', 'zoom', 'Units', 'normalized', 'Position', [0.91 0.32 0.0700 0.0333]);
            p.filterMasksThresh = uicontrol('Style', 'pushbutton', 'String', 'filter masked spots', 'Units', 'normalized', 'Position', [0.91 0.28 0.0700 0.0333]);
            p.threshValue = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [0.835 0.32 0.0700 0.0333]);

            p.upperContrastSlider = uicontrol('Style', 'slider', 'Value', 0.5, 'Units', 'normalized', 'Position', [0.8444 0.61 0.1233 0.0050]);
            p.lowerContrastSlider = uicontrol('Style', 'slider', 'Units', 'normalized', 'Position', [0.8444 0.58 0.1233 0.0050]);
            p.sliderLabel = uicontrol('Style', 'text', 'String', 'Adjust contrast', 'Units', 'normalized', 'Position', [0.8750 0.53 0.07 0.03]);

            p.figHandle.Visible = 'on';
        end
        
        %attach controller
        function startupFcn(p)
            %Do we want to reassign spots to cells and update all masks
            %when relaunching GUI? Helpful if tables aren't saved when shutting down. 
            %Otherwise, probably unnecessary. 
%             p.nucleiObj.updateAllMasks;
%             p.spotTable.assignSpotsToNuclei;
%             p.spotTable.updateAllMasks;
            
            p.mainAxesCntrlr = d2MainAxesController(p, p.scanObj, p.spotTable, p.maskObj, p.nucleiObj);  
            p.attatchMainAxesController(p.mainAxesCntrlr);        
            
            p.threshAxesCntrlr = d2ThresholdAxesController(p, p.scanObj, p.spotTable, p.maskObj, p.nucleiObj);
            p.attachThresholdController(p.threshAxesCntrlr);
             
            p.thumbAxesCntrlr = d2ThumbnailAxesController(p, p.scanObj, p.spotTable);  
            p.attachThumbnailController(p.thumbAxesCntrlr);
            
            %Connect controllers
            p.mainAxesCntrlr.threshCntrlr = p.threshAxesCntrlr;
            p.mainAxesCntrlr.thumbCntrlr = p.thumbAxesCntrlr;
            
            p.threshAxesCntrlr.mainAxesCntrlr = p.mainAxesCntrlr;
            p.threshAxesCntrlr.thumbCntrlr = p.thumbAxesCntrlr;
            p.threshAxesCntrlr.startup();
            
            p.thumbAxesCntrlr.mainAxesCntrlr = p.mainAxesCntrlr;
            p.thumbAxesCntrlr.threshCntrlr = p.threshAxesCntrlr;
            p.thumbAxesCntrlr.startup();
        end
        
        function p = attatchMainAxesController(p, controller)
            p.channelPopup.Callback = {@controller.changeChannel};
            p.colormapPopup.Callback = {@controller.changeColormap};
            p.centroidList.Callback = {@controller.centroidSelected};
            p.spotsCheckBox.Callback = {@controller.overlaySpots};
            p.centroidsCheckBox.Callback = {@controller.overlayNuclei};
            p.scatterCheckBox.Callback = {@controller.scatterCallback};
            p.masksCheckBox.Callback = {@controller.overlayMasks};
            p.addCellButton.Callback = {@controller.addCells};
            p.deleteCellButton.Callback = {@controller.deleteCells};
            p.maskSpotsButton.Callback = {@controller.addSpotMask};
            %p.maskSpotsAllChannelsCheckBox.Callback={@controller.toggleMaskSpotsAllChannels}
            p.maskCellButton.Callback = {@controller.addCellMask};
            p.deleteMaskButton.Callback = {@controller.deleteMask2};
            p.zoomAxes.Callback = {@controller.zoomInPressed};
            p.panAxes.Callback = {@controller.panViewPressed};
            p.upperContrastSlider.Callback = {@controller.updateMainAxes};
            p.lowerContrastSlider.Callback = {@controller.updateMainAxes};
            p.saveButton.Callback = {@p.saveButtonPressed};
            p.exportButton.Callback = {@p.exportButtonPressed};
            p.shuffleColors.Callback = {@controller.shuffleColorsInView};
            
            p.figHandle.WindowButtonDownFcn = {@controller.figWindowDown};
            p.figHandle.KeyPressFcn = {@controller.keyPressFunctions};
            p.figHandle.CloseRequestFcn = {@p.closeFigFcn}; 
            
            %Set KeyPressFcn for all uicontrols. Change this if you want
            %some uicontrols to have unique KeyPressFcn (e.g. arrow keys) 
            set(findobj(p.figHandle, 'Type', 'UIControl'), 'KeyPressFcn', {@controller.keyPressFunctions})
        end
        
        function p = attachThresholdController(p, controller)
            p.threshValue.Callback = {@controller.threshValueChange};
            p.threshAxes.ButtonDownFcn = {@controller.thresholdButtonDown};
            p.zoomThresh.Callback = {@controller.threshZoomButtonDown};
            p.filterMasksThresh.Callback = {@controller.filterMasksPushed};
            %p.threshValue.KeyPressFcn = {@controller.threshValueKey};
        end
        
        function p = attachThumbnailController(p, controller)
            p.thumbAxes.ButtonDownFcn = {@controller.thumbAxesButtonDown};
        end
        
        function saveButtonPressed(p, ~, ~)
            fprintf('Saving mask table, cell table, and spot tables.\nThis may take a minute\n')
            p.nucleiObj.saveNucleiTable;
            p.spotTable.updateScanSummary;
            p.spotTable.saveSpotsTable;
            p.maskObj.saveMasksTable;
            disp('done')
        end
        
        function exportButtonPressed(p, ~, ~)
            fprintf('Saving %s\n',p.spotTable.spotsSummaryFile)
            p.spotTable.exportSpotsSummary
            disp('done')
        end
        
        function relaunchGUI(p)
            createComponents(p)
            startupFcn(p)
        end
        
        function closeFigFcn(p, ~, ~)
            delete(p.figHandle)
            p.saveButtonPressed;
        end
        
    end

end
