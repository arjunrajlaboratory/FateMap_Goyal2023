%--------------------------------------------------------------------------
% This GUI displays a series of layouts and allows the user to select the
% appropriate one. Hs.layoutIndex is needed by getFilePaths  
%--------------------------------------------------------------------------
function Hs = getLayoutOrientation(Hs)
    Td = populateGUI;
    
    Td.rows = Hs.rows;
    Td.cols = Hs.cols;
    Td.nameExt = Hs.nameExt;
    Td.foundChannels = Hs.foundChannels;
    Td.dirPath = Hs.dirPath;   

    %Plot the grid
    set(Td.gridAx,'xlim',[0,3],'ylim',[0,3]);
    hold on;
    plot([0,3],[0,0]);
    plot([0,3],[1,1]);
    plot([0,3],[2,2]);
    plot([0,3],[3,3]);

    plot([0,0],[0,3]);
    plot([1,1],[0,3]);
    plot([2,2],[0,3]);
    plot([3,3],[0,3]);
    hold off;

    %Create array of cells containing the different layouts
    %Row no-snake
    lay1 = [1,2,3;4,5,6;7,8,9];
    %Row snake
    lay2 = [1,2,3;6,5,4;7,8,9];
    %Row-flipped no-snake
    lay3 = [3,2,1;6,5,4;9,8,7];
    %Row-flipped snake
    lay4 = [3,2,1;4,5,6;9,8,7];
    %Col no-snake
    lay5 = [1,4,7;2,5,8;3,6,9];
    %Col snake
    lay6 = [1,6,7;2,5,8;3,4,9];
    %Col-flipped no-snake
    lay7 = [3,6,9;2,5,8;1,4,7];
    %Col-flipped snake
    lay8 = [3,4,9;2,5,8;1,6,7];
    %Col&Row flipped no snake
    lay9 = [9,6,3;8,5,2;7,4,1];
    %Col&Row flipped snake
    lay10 = [9,4,3;8,5,2;7,6,1];

    Td.layouts = [mat2cell(lay1),mat2cell(lay2),mat2cell(lay3),mat2cell(lay4),...
        mat2cell(lay5),mat2cell(lay6),mat2cell(lay7),mat2cell(lay8),mat2cell(lay9),mat2cell(lay10)];
    Td.layoutStrings = [mat2cell('Row no snake'),mat2cell('Row snake'),mat2cell('Row flipped no snake'),...
        mat2cell('Row flipped snake'), mat2cell('Col no snake'), mat2cell('Col snake'),...
        mat2cell('Col flipped no snake'),mat2cell('Col flipped snake'),mat2cell('Col&Row flipped no snake'),mat2cell('Col&Row flipped snake')];
    Td = populateRadioButtons(Td);
    
    Td.layoutIndex = 1;
    Td = drawLayout(Td);
    guidata(Td.figH,Td);
    %waits untill uiresume called or until figure is deleted
    uiwait(Td.figH);
    %uiresume(Hs.figH) is called when select or exit called and then code
    %that is below is executed
    Td = guidata(Td.figH);
    Hs.layoutIndex = Td.layoutIndex;
    delete(Td.figH);
end
function Td = populateGUI
    figH = figure('Position',[400 400 520 300],...
        'NumberTitle','off',...
        'Name','Layout Selector',...
        'Resize','on',...
        'Toolbar','none',...
        'MenuBar','none',...
        'Color',[0.247 0.247 0.247],...
        'CloseRequestFcn',@closeRequestCallBack,...
        'Visible','on');

    Td = guihandles(figH);
    Td.figH = figH;
    Td.mainPanel = uipanel('Parent',Td.figH,...
        'Units','normalized',...
        'BorderType','etchedin',...
        'BackgroundColor',[0.247 0.247 0.247],...
        'Visible','on',...
        'Position',[0,0,0.577,1]);
    Td.rightPanel = uipanel('Parent',figH,...
        'Units','normalized',...
        'BorderType','etchedin',...
        'BackgroundColor',[0.247 0.247 0.247],...
        'Visible','on',...
        'Position',[0.587,0,0.413,1]);
    Td.buttonGroup = uibuttongroup('Parent',Td.rightPanel,...
        'Units','normalized',...
        'Visible','on',...
        'SelectionChangeFcn',@radioButtonCallBack,...
        'BackgroundColor',[0.247 0.247 0.247],...
        'Position',[0,0,1,1]);
    Td.gridAx = axes('Parent',Td.mainPanel,...
        'Units','normalized',...
        'Position',[0.05,0.3,.9,.67],...
        'YDir','reverse',...
        'XTick',[],'YTick',[],...
        'Color',[1,1,1]);
    Td.nextButton  = uicontrol('Parent',Td.mainPanel,...
        'String','Next',...
        'Style','pushbutton',...
        'FontSize',10,...
        'HorizontalAlignment','Left',...
        'Units','normalized',...
        'Position',[0.20 0.15 0.275 0.1],...
        'BackgroundColor',[1 1 1],...
        'Callback',@nextLayoutCallBack);
    Td.prevButton  = uicontrol('Parent',Td.mainPanel,...
        'String','Previous',...
        'Style','pushbutton',...
        'FontSize',10,...
        'HorizontalAlignment','Left',...
        'Units','normalized',...
        'Position',[0.20 0.05 0.275 0.1],...
        'BackgroundColor',[1 1 1],...
        'Callback',@prevLayoutCallBack);
    Td.previewButton  = uicontrol('Parent',Td.mainPanel,...
        'String','Preview',...
        'Style','pushbutton',...
        'FontSize',10,...
        'HorizontalAlignment','Left',...
        'Units','normalized',...
        'Position',[0.50 0.15 0.275 0.1],...
        'BackgroundColor',[1 1 1],...
        'Callback',@previewButtonCallBack);
    Td.selectButton  = uicontrol('Parent',Td.mainPanel,...
        'String','Select',...
        'Style','toggle',...
        'FontSize',10,...
        'HorizontalAlignment','Left',...
        'Units','normalized',...
        'Position',[0.50 0.05 0.275 0.1],...
        'BackgroundColor',[1 1 1],...
        'Callback',@selectAndCloseCallBack);
end
function previewButtonCallBack(hObject, eventData)
    Td = guidata(gcbo);
    Td = getFilePaths(Td);
    if isfield(Td,'figHPrev') && ishandle(Td.figHPrev)
        delete(Td.figHPrev);
    end
    Td.figHPrev = figure('Position',[500 200 700 650],...
        'NumberTitle','off',...
        'Name','Preview',...
        'Resize','on',...
        'Toolbar','none',...
        'MenuBar','none',...
        'Color',[0.247 0.247 0.247],...
        'CloseRequestFcn',@closeRequestCallBack,...
        'Visible','on');
    Td.bottomPanelPrev = uipanel('Parent',Td.figHPrev,...
        'Units','normalized',...
        'BorderType','etchedin',...
        'BackgroundColor',[0.247 0.247 0.247],...
        'Visible','on',...
        'Position',[0.05,0.01,0.86,0.10]);
    Td.imgAxPrev = axes('Parent',Td.figHPrev,...
        'Units','normalized',...
        'Position',[0.05,0.12,.86,.86],...
        'YDir','reverse',...
        'XTick',[],'YTick',[],...
        'Color',[1,1,1]);

    axis equal;

    Td.randomButtonPrev  = uicontrol('Parent',Td.bottomPanelPrev,...
        'String','Random',...
        'Style','pushbutton',...
        'FontSize',10,...
        'HorizontalAlignment','Left',...
        'Units','normalized',...
        'Position',[0.29 0.275 0.13 0.45],...
        'BackgroundColor',[1 1 1],...
        'Callback',@randomButtonPrevCallBack);
    Td.contrastButtonPrev  = uicontrol('Parent',Td.bottomPanelPrev,...
        'String','Contrast',...
        'Style','toggle',...
        'FontSize',10,...
        'HorizontalAlignment','Left',...
        'Units','normalized',...
        'Position',[0.59 0.275 0.13 0.45],...
        'BackgroundColor',[1 1 1],...
        'Callback',@contrastButtonPrevCallBack);
    Td.foundChannelsTemp = [mat2cell('dapi')];
    for channel = Td.foundChannels
        if ~strcmp(cell2mat(channel),'dapi')
            Td.foundChannelsTemp = [Td.foundChannelsTemp,channel];
        end
    end
    Td.chanPopPrev  = uicontrol('Parent',Td.bottomPanelPrev,...
        'Style','popup',...
        'String',Td.foundChannelsTemp,...
        'FontSize',10,...
        'HorizontalAlignment','Left',...
        'Units','normalized',...
        'Position',[0.44 0.275 0.13 0.45],...
        'ForegroundColor',[1 1 1],...
        'BackgroundColor',[0.247 0.247 0.247],...
        'Callback',@chanPopPrevCallBack);
    Td.currRowPrev = 1;
    Td.currColPrev = 1;
    Td.chanIndexPrev = 1;
    Td = displayImages(Td);
    guidata(Td.figH,Td);
    guidata(Td.figHPrev,Td);
end
function chanPopPrevCallBack(hObject,eventdata)
    Td = guidata(gcbo);
    Td = guidata(Td.figH);
    set(Td.contrastButtonPrev,'Value',0);
    Td.chanIndexPrev = get(hObject,'Value');
    Td = displayImages(Td);
    guidata(Td.figH,Td);
    guidata(Td.figHPrev,Td);
end
function contrastButtonPrevCallBack(hObject, eventData)
    Td = guidata(gcbo);
    Td = guidata(Td.figH);
    Td = displayImages(Td);
    guidata(Td.figH,Td);
    guidata(Td.figHPrev,Td);
end
function randomButtonPrevCallBack(hObject, eventData)
    % gcbo refers to Td.figHPrev which has a copy of Td.figH
    Td = guidata(gcbo);
    % However, the main data is stored with Td.figH
    Td = guidata(Td.figH);
    Td.currRowPrev = ceil(rand * (Td.rows - 1));
    Td.currColPrev = ceil(rand * (Td.cols - 1));
    Td = displayImages(Td);
    guidata(Td.figH,Td);
end
function Td = displayImages(Td)
    currImage = imread(cell2mat(Td.filePaths(Td.currRowPrev,Td.currColPrev,Td.chanIndexPrev)));
    Td.imageSize = size(currImage);
    minCurr = min(currImage(:));
    maxCurr = max(currImage(:));
    
    rightImage = imread(cell2mat(Td.filePaths(Td.currRowPrev,Td.currColPrev+1,Td.chanIndexPrev)));
    minRight = min(rightImage(:));
    maxRight = max(rightImage(:));
    
    downImage = imread(cell2mat(Td.filePaths(Td.currRowPrev+1,Td.currColPrev,Td.chanIndexPrev)));
    minDown = min(downImage(:));
    maxDown = max(downImage(:));
    
    downRightImage = imread(cell2mat(Td.filePaths(Td.currRowPrev+1,Td.currColPrev+1,Td.chanIndexPrev)));
    minDownRight = min(downRightImage(:));
    maxDownRight = max(downRightImage(:));
    
    imgCat = [currImage,rightImage;downImage,downRightImage];
    
    Td.imgCat = imgCat;
    
    set(Td.imgAxPrev,'xlim',[0,(Td.imageSize(2) * 2)],'ylim',[0,(Td.imageSize(1) * 2)]);

    minInt = median(double([minCurr,minDown,minRight,minDownRight]));
    maxInt = median(double([maxCurr,maxDown,maxRight,maxDownRight]));
    
    imgCatTemp = scale(imgCat,[minInt,maxInt]);
    
    if get(Td.contrastButtonPrev,'Value') == 1
        imgCatTemp = imgCatTemp * 2;
    end
    imshow(imgCatTemp,'Parent',Td.imgAxPrev);
end
function Td = populateRadioButtons(Td)
    Td.radioButtonsH = [];
    spacing = (1 - 0.02) / numel(Td.layouts);
    for index = 1:numel(Td.layouts)
        % yLoc starts at top and goes down.
        % (spacing - 0.01) is the height of a radio-button
        yLoc = spacing * numel(Td.layouts) - (spacing - 0.01) - (spacing * (index - 1));
        radioButton = uicontrol('Parent',Td.buttonGroup,...
            'Units','normalized',...
            'Style','radiobutton',...
            'BackgroundColor',[0.247 0.247 0.247],...
            'FontSize',10,...
            'Visible','on',...
            'ForegroundColor',[1,1,1],...
            'String',cell2mat(Td.layoutStrings(index)),...
            'Position',[0,yLoc,0.9,spacing - 0.01]); 
        Td.radioButtonsH = [Td.radioButtonsH,radioButton];
    end
end
function radioButtonCallBack(hObject, eventData)
    Td = guidata(gcbo);
    str = get(eventData.NewValue,'String');
    val = 1;
    for index = 1:numel(Td.layoutStrings)
        layout =  cell2mat(Td.layoutStrings(index));
        if strcmp(str,layout)
            val = index;
            break;
        end
    end
    Td.layoutIndex = val;

    Td = drawLayout(Td);
    guidata(Td.figH,Td);    
end
function selectAndCloseCallBack(hObject, events)
    Td = guidata(gcbo);
    uiresume(Td.figH);
end
function Td = drawLayout(Td)
    if isfield(Td,'text')
        for p = Td.text
            if ishandle(p)
                delete(p);
            end
        end
    end
    layout = cell2mat(Td.layouts(Td.layoutIndex));
    rLoc = 0.5;
    cLoc = 0.5;
    hold on;
    rLoc = 0.5;
    Td.text = [];
    for row = 1:3
        cLoc = 0.5;
        for col = 1:3
            Td.text = [Td.text,text(cLoc,rLoc,int2str(layout(row,col)))];
            cLoc = cLoc + 1;
        end
        rLoc = rLoc + 1;
    end    
end
function closeRequestCallBack(hObject,events)
    if isequal(get(hObject,'waitstatus'),'waiting')
        uiresume(hObject);
    else
        delete(hObject);
    end
end
function nextLayoutCallBack(hObject,events)
    Td = guidata(gcbo);
    Td.layoutIndex = Td.layoutIndex + 1;  
    if Td.layoutIndex > size(Td.layouts)
        Td.layoutIndex = 1;
    end
    set(Td.radioButtonsH(Td.layoutIndex),'Value',1);
    Td = drawLayout(Td);
    guidata(Td.figH,Td);
end
function prevLayoutCallBack(hObject,events)
    Td = guidata(gcbo);
    Td.layoutIndex = Td.layoutIndex - 1;  
    if Td.layoutIndex > size(Td.layouts)
        Td.layoutIndex = 1;
    elseif Td.layoutIndex < 1
        Td.layoutIndex = numel(Td.layouts);
    end
    set(Td.radioButtonsH(Td.layoutIndex),'Value',1);

    Td = drawLayout(Td);
    guidata(Td.figH,Td);
end    