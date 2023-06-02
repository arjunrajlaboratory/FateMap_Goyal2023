classdef d2ThresholdGUI < handle
    
    properties (Access = public)
        
        spotTable
        scanObject
        
        addCurrButtonHandle
        addNextButtonHandle
        deleteButtonHandle
        deselectAllButtonHandle
        connectParentButtonHandle
        saveButtonHandle
        showPointIDsHandle
        currentFramePopupHandle
        
        noneButtonHandle
        badButtonHandle
        deadButtonHandle
        disappearButtonHandle
        appearButtonHandle
        missingButtonHandle
        leftFrameButtonHandle
        enterFrameButtonHandle
        
        axesHandle
        imageHandle
        
    end
    
    methods
        function p = d2ThresholdGUI(spotTable,scanObject)
            p.figHandle = figure('Visible','off','Position',[360,500,450,285]);
            %p.figHandle.KeyPressFcn = @p.GUIWindowKeyPressFcn;
            
            p.addCurrButtonHandle = uicontrol('Style','pushbutton','String','add curr pt','Position',[5,220,80,25]);
            p.addNextButtonHandle = uicontrol('Style','pushbutton','String','add next pt','Position',[5,250,80,25]);
            p.deleteButtonHandle = uicontrol('Style','pushbutton','String','delete pts','Position',[5,280,80,25]);
            p.deselectAllButtonHandle = uicontrol('Style','pushbutton','String','deselect','Position',[5,310,80,25]);
            p.connectParentButtonHandle = uicontrol('Style','pushbutton','String','connect parent','Position',[5,340,80,25]);
            p.saveButtonHandle = uicontrol('Style','pushbutton','String','save','Position',[5,370,80,25]);
            p.showPointIDsHandle = uicontrol('Style','checkbox','String','showID','Position',[5,400,80,25]);
            
            p.noneButtonHandle =        uicontrol('Style','pushbutton','String','none','Position',[5,430,80,25]);
            p.badButtonHandle =         uicontrol('Style','pushbutton','String','bad','Position',[5,460,80,25]);
            p.deadButtonHandle =        uicontrol('Style','pushbutton','String','dead','Position',[5,490,80,25]);
            p.disappearButtonHandle =   uicontrol('Style','pushbutton','String','disappear','Position',[5,520,80,25]);
            p.appearButtonHandle =      uicontrol('Style','pushbutton','String','appear','Position',[5,550,80,25]);            
            p.missingButtonHandle =     uicontrol('Style','pushbutton','String','missing','Position',[5,580,80,25]);
            p.leftFrameButtonHandle =   uicontrol('Style','pushbutton','String','leftFrame','Position',[5,610,80,25]);
            p.enterFrameButtonHandle =  uicontrol('Style','pushbutton','String','enterFrame','Position',[5,640,80,25]);
            
            p.currentFramePopupHandle = uicontrol(p.figHandle,'Style','popupmenu');
            p.currentFramePopupHandle.Position = [5 5 200 25];
            
            
            
            p.axesHandle = axes('Units','normalized','Position',[0.05 0.02 .98 .98]);
            
            imshow(rand(100),'Parent',p.axesHandle);
            
            % Maybe keep invisible until we load data?
            p.figHandle.Visible = 'on';
            
            
            % Should fix this with a proper display method
            
            % p.pointController = pointController(p.pointTableHandle, p.addCurrButtonHandle, p.addNextButtonHandle, p.deselectAllButtonHandle, p.deleteButtonHandle, p.connectParentButtonHandle, p.currentFramePopupHandle, p.axesHandle, p.linesHandle);
            p.pointController = pointController(p.pointTableHandle, p.addCurrButtonHandle, p.addNextButtonHandle, p.deselectAllButtonHandle, p.deleteButtonHandle, p.connectParentButtonHandle, p.currentFramePopupHandle, p.axesHandle, p.linesHandle);
            pc = p.pointController;

            p.addCurrButtonHandle.Callback={@pc.addCurrPointButtonPushed}; % Needs fixing
            p.addNextButtonHandle.Callback={@pc.addNextPointButtonPushed}; % Needs fixing
            p.deleteButtonHandle.Callback ={@pc.deleteButtonPushed}; % Needs fixing
            p.deselectAllButtonHandle.Callback ={@pc.deselectAllButtonPushed}; % Needs fixing
            p.connectParentButtonHandle.Callback ={@pc.connectParentButtonPushed}; % Needs fixing
            p.currentFramePopupHandle.Callback ={@pc.updateFrame}; % Needs fixing
            p.saveButtonHandle.Callback ={@pc.saveButtonPushed}; % Needs fixing
            p.showPointIDsHandle.Callback ={@pc.showPointIDsPushed}; % Needs fixing
            
            p.noneButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.badButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.deadButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.disappearButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.appearButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.missingButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.leftFrameButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            p.enterFrameButtonHandle.Callback ={@pc.setAnnotation}; % Needs fixing
            
            pc.showPointIDsHandle = p.showPointIDsHandle;
            pc.saveFilename = inFilename;
            

            p.loadAllData();
            pc.pointTableHandle = p.pointTableHandle;
            p.pointController.saveButtonPushed([],[]);

            %pc.makePoints();
            pc.showImages();
            
        end
        
        
        function p = loadAllData(p) % Lots of this should probably be in a fileController
            
            %%% SHOULD USE THIS
            p.fileTable = parseFiles();
            %Probably should call writetable within parseFiles
            writetable(p.fileTable,'fileTable.csv');
            
            % Keep only the histone marker channel
            fileTable = p.fileTable(p.fileTable.wavelength == "cy5",:);
            
            fileTable = sortrows(fileTable,'frameNumber');
            
            % Replace with some sort of file loader
            %p.imageFiles = {'Scan210_w1_s18_t1.TIF','Scan211_w1_s18_t1.TIF','Scan210_w1_s18_t1.TIF'};
            p.imageFiles = fileTable.fileName;
            
            
            p.currentFramePopupHandle.Value = 1;
            p.currentFramePopupHandle.String = p.imageFiles(1:end-1);
            p.currentFramePopupHandle.UserData = p.imageFiles;
            
            if isfile(p.pointController.saveFilename)
                p.pointTableHandle = pointTable(p.pointController.saveFilename);
                %p.pointTableHandle.allPoints = readtable(p.pointController.saveFilename,'TextType','string');
            else
                
                p.pointTableHandle = pointTable();
                allCents = zeros(0,2);
                allIntensities = zeros(0);
                allFrames = zeros(0);
                for i = 1:length(p.imageFiles)
                    fprintf('Reading and finding cells in image %d\n',i);
                    currImage = imread(p.imageFiles{i});
                    %cents = getCentroidsTest(currImage);
                    %p.pointTableHandle.addRawPoints(i,cents);
                    
                    
                    [cents, intensities] = getCentroids(currImage);
                    allCents = [allCents; cents];
                    allIntensities = [allIntensities intensities];
                    allFrames = [allFrames i*ones(1,length(intensities))];
                end
                
                % create a histogram of the centroid intensities:
                figure('visible', 'off');
                histogram_handle = histogram(allIntensities);
                histogram_x = histogram_handle.BinEdges;
                histogram_y = histogram_handle.Values;
                
                % get the inflection point of the histogram:
                histogram_y_derivative = diff(histogram_y);
                inflection_y = find(histogram_y_derivative>0, 1) + 1;
                inflection_x = 1/2 * (histogram_x(inflection_y) + histogram_x(inflection_y + 1));
                
                % keep the centroids past the inflection point of the histogram:
                idx = allIntensities > inflection_x;
                
                allCents = allCents(idx,:);
                allFrames = allFrames(idx);
                
                % Add in all the points to the pointTable
                for i = 1:length(p.imageFiles)
                    idx = allFrames == i;
                    cents = allCents(idx,:);
                    p.pointTableHandle.addRawPoints(i,cents);
                end
                
                
                p.pointTableHandle.guessParents();
            end
            
            
            % Add in something to draw the cells. Actually, no, let the
            % earlier call deal with this.
            
            %pc = p.pointController;
            %pc.showImages();

            
            
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
end
