classdef d2ThresholdView < matlab.apps.AppBase
    
    properties (Access = public)
        
        %Model objects
        spotTable
        scanObj
        maskObj
        nucleiObj
        
        controlObj
        
        figHandle
        mainAxes 
        thumbAxes
        threshAxes
        centroidList %List box 
        spotsCheckBox %Check box
        centroidsCheckBox %Check box
        

        channelPopup %Drop down
        saveButton %Push button
        exportButton %Push button
        homeButton %Push button
        
        zoomInThresh %Push button
        zoomOutThresh %Push button
        resetThresh %Push button
        
        addCellButton %Toggle button
        deleteCellButton %Toggle button
        maskCellButton %Push button
        maskSpotsButton %Push button
        deleteMaskButton %Push button
        
        
        upperContrastSlider %Slider
        lowerContrastSlider %Slider
        sliderLabel
    end
    
    methods (Access = private)
        
        % Create UIFigure and components
        function createComponents(app)

            app.figHandle = uifigure('Visible', 'off', 'Position', [200 100 900 605]);

            app.mainAxes = uiaxes(app.figHandle, 'XTickLabel', '', 'YTickLabel', '', 'Position', [1 1 550 550]);
            app.thumbAxes = uiaxes(app.figHandle, 'XTickLabel', '', 'YTickLabel', '', 'Position', [549 1 331 210]);
            app.threshAxes = uiaxes(app.figHandle, 'XTickLabel', '', 'YTickLabel', '', 'Position', [549 361 190 190]);
            
            app.centroidList = uilistbox(app.figHandle, 'Position', [758 381 115 151]);
            app.DropDown = uidropdown(app.figHandle, 'Items', app.spotTable.spotChannels, 'Position', [565 304 100 22]);          
            app.spotsCheckBox = uicheckbox(app.figHandle, 'Text', 'spots', 'Position', [573 347 65 20]);
            app.centroidsCheckBox = uicheckbox(app.figHandle, 'Text', 'nuclei', 'Position', [651 347 65 20]);

            app.addCellButton = uibutton(app.figHandle, 'state', 'Text', 'add cells', 'Position', [781 242 100 22]);
            app.deleteCellButton = uibutton(app.figHandle, 'state', 'Text', 'delete cells', 'Position', [781 210 100 22]);

            app.maskSpotsButton = uibutton(app.figHandle, 'push', 'Text', 'mask spots', 'Position', [565 242 100 22]);
            app.maskCellButton = uibutton(app.figHandle, 'push', 'Text', 'mask cells', 'Position', [671 210 100 22]);
            app.deleteMaskButton = uibutton(app.figHandle, 'push', 'Text', 'delete mask', 'Position', [671 242 100 22]);
            
            app.homeButton = uibutton(app.figHandle, 'push', 'Text', 'home', 'Position', [565 210 100 22]);
            app.saveButton = uibutton(app.figHandle, 'push', 'Text', 'save', 'Position', [684 304 70 22]);
            app.exportButton = uibutton(app.figHandle, 'push', 'Text', 'export', 'Position', [684 274 70 22]);
            
            app.zoomInThresh = uibutton(app.figHandle, 'push', 'Text', 'zoom in', 'Position', [806 145 63 20]);
            app.zoomOutThresh = uibutton(app.figHandle, 'push', 'Text', 'zoom out', 'Position', [806 122 63 20]);
            app.resetThresh = uibutton(app.figHandle, 'push', 'Text', 'reset', 'Position', [806 169 63 20]);

            app.upperContrastSlider = uislider(app.figHandle, 'MajorTicks', [], 'MajorTickLabels', {''}, 'MinorTicks', [], 'Position', [760 355 111 3]);
            app.lowerContrastSlider = uislider(app.figHandle, 'MajorTicks', [], 'MajorTickLabels', {''}, 'MinorTicks', [], 'Position', [760 330 111 3]);
            app.sliderLabel = uilabel(app.figHandle, 'Text', 'Adjust contrast', 'Position', [774 299 87 22]);

            app.UIFigure.Visible = 'on';
        end
        
        %attach controller
        function startupFcn(app)
            app.controlObj = Controller(app, app.scanObj, app.spotTable, app.maskObj, app.nucleiObj);  
            app.attatchToController(app.controlObj);        
            
            %app.modelObj.addlistener('balanceChanged',@app.updateBalance);  
        end
        
        function attatchToController(app, controller)
            %funcH = @controller.callback_withDrawButton;
            addlistener(app.maskSpotsButton,'ButtonPushed',funcH)

            addlistener(app.maskCellButton,'ButtonPushed',funcH) 
            
            addlistener(app.deleteMaskButton,'ButtonPushed',funcH)
            
            addlistener(app.homeButton,'ButtonPushed',funcH) 
            
            addlistener(app.saveButton,'ButtonPushed',funcH)
            
            addlistener(app.exportButton,'ButtonPushed',funcH) 
            
            addlistener(app.addCellButton,'ValueChanged',funcH)
            
            addlistener(app.centroidList,'ValueChanged',funcH) 
            
            addlistener(app.DropDown,'ValueChanged',funcH) 
            
            addlistener(app.upperContrastSlider,'ValueChanged',funcH) 
            
            addlistener(app.lowerContrastSlider,'ValueChanged',funcH)
            
            addlistener(app.zoomInThresh,'ButtonPushed',funcH) 
            
            addlistener(app.zoomOutThresh,'ButtonPushed',funcH)
            
            addlistener(app.resetThresh,'ButtonPushed',funcH) 
        end

    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = d2ThresholdView(scanObject, spotTable, maskObject)
            app.scanObj = scanObject;
            app.spotTable = spotTable;
            app.maskObj = maskObject;

            % Create UIFigure and components
            createComponents(app)
            
            registerApp(app, app.figHandle)
            startupFcn(app)
            %runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)
            %save tables

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
    
    % Component initialization
    methods (Access = private)
        function createComponents(p)
            p.figHandle = uifigure('Visible','off','Position',[200 100 900 605]);
            p.mainAxes = uiaxes(p.figHandle, 'Position', [1 1 550 550], 'XTickLabel', '', 'YTickLabel', '');
            p.thumbAxes = uiaxes(p.figHandle, 'Position', [549 1 331 210], 'XTickLabel', '', 'YTickLabel', '');
            p.threshAxes = uiaxes(p.figHandle, 'Position', [549 361 190 190], 'XTickLabel', '', 'YTickLabel', '');
            
            p.centroidListHandle = uilistbox(p.figHandle, 'Position', [573 347 67 22]);
            
            
            p.figHandle.Visible = 'on';
        end
    end
    
    methods (Access = public)
        function p = d2ThresholdGUI();
            p.scanObj = scanObject;
            p.spotTable = spotTable;
            p.maskObj = maskObject;
            createComponents(p)
            registerApp(p, p.figHandle);
            
        end
        p = d2ThresholdGUI(scanObject, spotTable, maskObject);
        createComponents(p
    end
    

        
        function p = loadData(p) % Need to make this a proper loader/analyzer
            
            image1 = imread('Scan210_w1_s18_t1.TIF');
            image2 = imread('Scan211_w1_s18_t1.TIF');
            
            cents1 = getCentroidsTest(image1);
            cents2 = getCentroidsTest(image2);
            
            outRGB = makeColoredImage(scale(im2double(image1)),[0 0.6797 0.9336]) + makeColoredImage(scale(im2double(image2)),[0.9648 0.5781 0.1172]);

            p.imageHandle = imshow(outRGB,'Parent',p.axesHandle);
            
            p.pointTableHandle = pointTable();
            
            p.pointTableHandle.addRawPoints(1,cents1);
            p.pointTableHandle.addRawPoints(2,cents2);
            
            p.pointTableHandle.guessParents();
            
        end
        
        
        function p = GUIWindowKeyPressFcn(p, src, eventdata)
            % determine the key that was pressed
            keyPressed = eventdata.Key;
            switch(keyPressed)
                case 'a'
                    p.pointController.addNextPointButtonPushed(src, eventdata);
                case 's'
                    p.pointController.addCurrPointButtonPushed(src, eventdata);
                case 'd'
                    p.pointController.deselectAllButtonPushed(src, eventdata);
                case 'f'
                    p.pointController.deleteButtonPushed(src, eventdata);
                case 'g'
                    p.pointController.connectParentButtonPushed(src, eventdata);
                case 'c'
                    p.currentFramePopupHandle.Value = max([1 p.currentFramePopupHandle.Value-1]);
                    p.pointController.updateFrame(src, eventdata);
                case 'v'
                    p.currentFramePopupHandle.Value = min([length(p.currentFramePopupHandle.String) p.currentFramePopupHandle.Value+1]);
                    p.pointController.updateFrame(src, eventdata);
                case 'z'
                    p.pointController.zoomMode();
                case 'q'
                    p.pointController.zoomReset();                    
                case 'x'
                    p.pointController.unZoom();
                case 't'
                    p.pointController.showPointIDsHandle.Value = ~p.pointController.showPointIDsHandle.Value;
                    p.pointController.showPointIDsPushed(p.pointController.showPointIDsHandle,eventdata);
                case 'w'
                    p.pointController.toggleGFP(src,eventdata);
                case 'e'
                    p.pointController.toggleTrans(src,eventdata);
            end
            
        end
end
