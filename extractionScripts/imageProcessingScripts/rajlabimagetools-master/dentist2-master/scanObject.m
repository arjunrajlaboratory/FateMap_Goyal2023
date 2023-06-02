classdef scanObject < handle
    
    properties (Access = public)
        
        scanFile
        scanSummaryFile = 'scanSummary.txt'
        tilesTableName = 'tilesTable.csv'

        tilesTable
        tileSize
        channels
        channelTypes
        scanMatrix
        scanDim
        snake
        startPos
        direction
        stitchDim
        dapiStitch
        smallDapiStitch
        stitchedScans
        smallStitchedScans
        resizeFactor = 4
        rowTransformCoords
        columnTransformCoords
        backgroundMats
        %Including imRotation in case the image acquired by Elements is flipped or somehow
        %rotated. Consider deleting. Also, height and width 
        %may need to be switched for rotated images. I'm not sure if 
        %tilesize changes when the image from Elements is rotated.
        imRotation = 0;
        prestitchedScanFileList
    end
    
    methods
        
        function p = scanObject(varargin) % 
            n = inputParser;
            n.addParameter('scanSummary', '', @ischar);
            n.addParameter('scanFile', '', @ischar); 
            n.addParameter('prestitchedScanFileList', {}, @(x) validateattributes(x,{'cell'},{'size',[1 nan]})); %Probably more accurate to have isempty(x) || (isvector(x) && iscellstr(x))
            n.addParameter('scanDim', [0,0], @(x)validateattributes(x,{'numeric'},{'size',[1 2]}));
            n.addParameter('channels',{}, @(x) isempty(x) || (isvector(x) && iscellstr(x))) 
            n.addParameter('channelTypes',{},@(x) all(ismember(x, {'dapi', 'FISH', 'other'})))
            n.parse(varargin{:});
            
            if ~isempty(n.Results.scanSummary) % a scan summary (Eg. scanSummary.txt) was given, load it
                    fprintf('Loading %s\n', n.Results.scanSummary);
                    p.scanSummaryFile = n.Results.scanSummary;
                    p.loadScanSummary(n.Results.channels);
                    p.loadChannelTypes(n.Results.channelTypes) % user-input would override the channelTypes from scanSummaryFile
            else
                if ~isempty(n.Results.scanFile) && all(n.Results.scanDim)
                    fprintf('Loading file %s\n', n.Results.scanFile)
                    p.scanFile = n.Results.scanFile;
                    p.channels = d2utils.loadChannelsND2(p.scanFile, n.Results.channels);
                    p.loadChannelTypes(n.Results.channelTypes)
                    p.scanDim = n.Results.scanDim;
                    p.defaultScanParameters()
                elseif ~isempty(n.Results.scanFile) %If loading pre-stitched scan. 
                    fprintf('Loading pre-stiched scan file: %s\n', n.Results.scanFile)
                    p.scanFile = n.Results.scanFile;
                    [~,~,ext] = fileparts(p.scanFile);
                    if strcmpi(ext, '.nd2')
                        p.channels = d2utils.loadChannelsND2(p.scanFile, n.Results.channels);
                        p.loadChannelTypes(n.Results.channelTypes)
                        p.loadPrestitchedScans();
                    elseif any(strcmpi(ext, {'.tif', '.tiff'}))
                        p.channels = d2utils.loadChannelsTiff(p.scanFile, n.Results.channels);
                        p.loadChannelTypes(n.Results.channelTypes)
                        p.loadPrestitchedScansTiff();
                    end
                elseif ~isempty(n.Results.prestitchedScanFileList)
                    temp=join(n.Results.prestitchedScanFileList,'  ');
                    fprintf('Loading pre-stitched scan file list:%s\n',temp{1})
                    p.scanFile='';
                    p.prestitchedScanFileList=n.Results.prestitchedScanFileList;
                    if isempty(n.Results.channels)
                        p.channels=replace(p.prestitchedScanFileList,'.tif','');
                    else
                        p.channels = n.Results.channels;
                    end
                    p.loadChannelTypes(n.Results.channelTypes)
                    p.loadPrestitchedScansFromFilelist()
                else
                    fprintf('Unable to create the scan object.\nPlease specify a scan summary file (e.g. scanSummary.txt) or the scan file name and scan dimensions.\n')
                    return
                end
            end
            
        end
        
        function defaultScanParameters(p)
            p.getTileSize();
            p.snake = true;
            p.direction = 'horizontal';
            p.startPos = 'top left';
            p.scanMatrix = d2utils.makeScanMatrix(p.scanDim);
        end
        
        function p = loadTiles(p)  
            
            height = p.tileSize(1);
            width = p.tileSize(2);
            nTiles = numel(p.scanMatrix);
            topCoords = zeros(nTiles,1);
            leftCoords = zeros(nTiles,1);
            for i = 1:nTiles
                [row,col] = find(p.scanMatrix == i);
                topCoords(i)  = col*p.columnTransformCoords(1) + row*p.rowTransformCoords(1);
                leftCoords(i) = row*p.rowTransformCoords(2) + col*p.columnTransformCoords(2);
            end
            
            topCoords = topCoords - min(topCoords) + 1;
            leftCoords = leftCoords - min(leftCoords) + 1;
            p.tilesTable = table((1:nTiles)', topCoords, leftCoords, repmat(height, nTiles,1), repmat(width, nTiles,1), ...
                'VariableNames', {'tileID', 'top', 'left', 'height', 'width'});
            p.stitchDim = [max(leftCoords)+p.tileSize(1)-1, max(topCoords)+p.tileSize(2)-1];
        end
          
         function p = getTileSize(p)
             reader = bfGetReader(p.scanFile);
             omeMeta = reader.getMetadataStore();
             p.tileSize = [omeMeta.getPixelsSizeY(0).getValue(), omeMeta.getPixelsSizeX(0).getValue()];
         end
         
         function outIm = getTileFromScan(p, tile, channel)
             reader = bfGetReader(p.scanFile);
             reader.setSeries(tile-1)
             channelIdx = find(ismember(p.channels, channel));  
             iPlane = reader.getIndex(0, channelIdx-1, 0) + 1;
             outIm  = bfGetPlane(reader, iPlane);
         end
         
         function p = loadPrestitchedScans(p)
             reader = bfGetReader(p.scanFile);
             reader.setSeries(0); %Will only load the first scan 
             %Try loading dapi. Assumes there is 1 dapi channel to load
             channelIdx = find(ismember(p.channelTypes, 'dapi'), 1);
             if ~isempty(channelIdx)
                 iPlane = reader.getIndex(0, channelIdx - 1, 0) + 1;
                 fprintf('   bfGetPlane: dapi\n')
                 p.dapiStitch = bfGetPlane(reader, iPlane);
             end
             %load other channels
             nonDapiChannels=p.channels(~ismember(p.channelTypes, 'dapi'));
             if ~isempty(nonDapiChannels)
                 p.stitchedScans.labels = nonDapiChannels;
                 p.stitchedScans.stitches = cell(1,numel(nonDapiChannels));
                 for i = 1:numel(nonDapiChannels)
                     channelIdx = find(ismember(p.channels, nonDapiChannels{i}));
                     iPlane = reader.getIndex(0, channelIdx - 1, 0) + 1;
                     fprintf('   bfGetPlane: %s (%i of %i non-dapi channels) with channelType=%s\n',nonDapiChannels{i},i,numel(nonDapiChannels),p.channelTypes{channelIdx})
                     p.stitchedScans.stitches{i} = bfGetPlane(reader, iPlane);
                 end
             end
             reader.close()
             p.scanDim = [1 1];
             p.stitchDim = size(p.dapiStitch);
         end
         
         function p = loadPrestitchedScansTiff(p)
             %first load dapi. Assumes there is 1 dapi channel to load
             dapiIdx = find(strcmpi('dapi', p.channelTypes), 1);
             if ~isempty(dapiIdx)
                 p.dapiStitch = imread(p.scanFile, 'Index', dapiIdx);
             end
             %load other channels
             nonDapiChannels=p.channels(~strcmpi('dapi', p.channelTypes));
             if ~isempty(nonDapiChannels)
                 p.stitchedScans.labels = nonDapiChannels;
                 p.stitchedScans.stitches = cell(1,numel(nonDapiChannels));
                 for i = 1:numel(nonDapiChannels)
                     channelIdx = find(strcmpi(nonDapiChannels{i}, p.channels));
                     fprintf('   bfGetPlane: %s (%i of %i non-dapi channels) with channelType=%s\n',nonDapiChannels{i},i,numel(nonDapiChannels),p.channelTypes{channelIdx})
                     p.stitchedScans.stitches{i} = imread(p.scanFile, 'Index', channelIdx);
                 end
             end
             p.scanDim = [1 1];
             p.stitchDim = size(p.dapiStitch);
         end
         
         function p = loadPrestitchedScansFromFilelist(p) % Ian D added
             
             %first load dapi
             dapiFile=p.prestitchedScanFileList{ismember(p.channelTypes,'dapi')};
             fprintf('   loading dapi from file %s\n',dapiFile)
             p.dapiStitch=imread(dapiFile);
             %load other channels
             nonDapiFilenames=p.prestitchedScanFileList(~ismember(p.channelTypes,'dapi'));
             nonDapiChannels=replace(nonDapiFilenames,'.tif','');
             p.stitchedScans.labels = nonDapiChannels;
             p.stitchedScans.stitches = cell(1,numel(nonDapiChannels));
             for i = 1:numel(nonDapiChannels)
                 fprintf('   loading %s (%i of %i non-dapi channels) from file %s\n',nonDapiChannels{i},i,numel(nonDapiChannels),nonDapiFilenames{i})
                 p.stitchedScans.stitches{i}=imread(nonDapiFilenames{i});
             end
             p.scanDim = [1 1];
             p.stitchDim = size(p.dapiStitch);
         end
        
         function loadChannelTypes(p,userInputChannelTypes) % Ian D added
             if isempty(userInputChannelTypes)
             % make default channelTypes based on filenames.
             %      anything with 'dapi' (case insensitive) in name = dapi
             %      'Brightfield' or 'trans' (case insensitive) = 'other'
             %      otherwise it is 'FISH'
                indOther=[];
                indDapi=find(contains(lower(p.channels),'dapi'));
                if length(indDapi)>1 && sum(ismember(p.channels,'dapi'))==1 % but there's just one exact match of 'dapi', use this and channelType of the rest 'other'
                    indOther=find(contains(lower(p.channels),'dapi') & ~ismember(p.channels,'dapi'));
                    indDapi=find(ismember(p.channels,'dapi'));
                else
                    if length(indDapi)>1
                        fprintf('WARNING: more than one channel was detected with dapi in the channel name (case insensitive). The first one (%s) will be used as the dapi channel, and the others will be assigned channelType other\n',p.channels{indDapi(1)})
                        indOther=indDapi(2:end);
                        indDapi=indDapi(1);
%                     elseif isempty(indDapi) %Ben updated 10/20/21 - allow for scans with no 'dapi' channels 
%                         error("no channel with dapi in the name was found. Try to input a channelTypes cell array, with a type for each channel. For a 5 channel file, it would be like this: launchD2ThresholdGUI(__,'channelTypes',{'dapi','FISH','FISH','other','FISH'}")
                    end
                end
                indTrans=find(contains(lower(p.channels),{'trans','brightfield'}));
                p.channelTypes=repmat({'FISH'},1,length(p.channels));
                p.channelTypes(indDapi)={'dapi'};
                p.channelTypes(indTrans)={'other'};
                p.channelTypes(indOther)={'other'}; % where multiple files with dapi (case insensitive) were found
             else % have user-input channelTypes
                p.checkOverChannelTypes(userInputChannelTypes)
                p.channelTypes=userInputChannelTypes;
             end
%              if isempty(userInputChannelTypes)
%                  if isfile(p.scanSummaryFile) 
%                      scanSummaryTable = d2utils.parseScanSummary(p.scanSummaryFile);
%                      p.channelTypes = split(scanSummaryTable{'channelTypes',1})';
%                      p.checkOverChannelTypes(p.channelTypes)
%                      
%                  elseif isempty(p.channelTypes)
%                      % make default channelTypes based on filenames.
%                      %      anything with 'dapi' (case insensitive) in name = dapi
%                      %      'Brightfield' or 'trans' (case insensitive) = 'other'
%                      %      otherwise it is 'FISH'
%                      indOther=[];
%                      indDapi=find(contains(lower(p.channels),'dapi'));
%                      if length(indDapi)>1 && sum(ismember(p.channels,'dapi'))==1 % but there's just one exact match of 'dapi', use this and channelType of the ris 'other'
%                          indOther=find(contains(lower(p.channels),'dapi') & ~ismember(p.channels,'dapi'));
%                          indDapi=find(ismember(p.channels,'dapi'));
%                      else
%                          if length(indDapi)>1
%                              fprintf('WARNING: more than one channel was detected with dapi in the channel name (case insensitive). The first one (%s) will be used as the dapi channel, and the others will be assigned channelType other\n',p.channels{indDapi(1)})
%                              indOther=indDapi(2:end);
%                              indDapi=indDapi(1);
%                          elseif isempty(indDapi)
%                              error("no channel with dapi in the name was found. Try to input a channelTypes cell array, with a cell for each channel. For a 5 channel file, it would be like this: launchD2ThresholdGUI(__,'channelTypes',{'dapi','FISH','FISH','other','FISH'}")
%                          end
%                      end
%                      indTrans=find(contains(lower(p.channels),{'trans','brightfield'}));
%                      p.channelTypes=repmat({'FISH'},1,length(p.channels));
%                      p.channelTypes(indDapi)={'dapi'};
%                      p.channelTypes(indTrans)={'other'};
%                      p.channelTypes(indOther)={'other'}; % where multiple files with dapi (case insensitive) were found
%                   % else p.channelTypes already exists, which is
%                   % unexpected.
%                  end
%              else % have user-input channelTypes
%                  p.checkOverChannelTypes(userInputChannelTypes)
%                  p.channelTypes=userInputChannelTypes;
%              end
         end
         
         function checkOverChannelTypes(p,channelTypes) % Ian D added
             mustBeMember(channelTypes,{'dapi','FISH','other'})
             if length(channelTypes)~=length(p.channels)
                 error('There are %i cells in channelTypes but there are %i channels',length(channelTypes),length(p.channels))
             end
             if sum(ismember(channelTypes,'dapi'))>1
                 error("There should be 1 dapi in the input channelTypes, but instead you have %i. If you actually do have more than one dapi channel, you must pick one to be processed as dapi and the rest can have channelType 'other'",sum(ismember(channelTypes,'dapi')))
%              elseif sum(ismember(channelTypes,'dapi'))==0 %Ben commented
%              out to allow loading scans without nuclei for d2IF, for example.
%                  error("One of the elements of channelTypes should be 'dapi'")
             end
         end
         
        %Stitch DAPI
        function p = stitchDAPI(p)  
            %Currently assumes there is only 1 dapi channel to stitch
            tileTable = p.tilesTable;
            tilesTmp = transpose(p.scanMatrix);
            tiles = tilesTmp(:);
            height = p.tileSize(1);
            width = p.tileSize(2);
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
            channel = find(ismember(p.channels, 'dapi')); 
%             channel = p.channels(ismember(p.channelTypes,'dapi'));
            reader = bfGetReader(p.scanFile);
            iPlane = reader.getIndex(0, channel - 1, 0) + 1;
            
            if logical(p.imRotation) %In case the image acquired by Elements is flipped or somehow rotated. Consider deleting. 
                for i = 1:numel(tiles)
                    reader.setSeries(tiles(i)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);
                    tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                        tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = rot90(tmpPlane, p.imRotation);
                end
            else
                for i = 1:numel(tiles)
                    reader.setSeries(tiles(i)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);
                    tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                        tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = tmpPlane;
                end
            end
            
            reader.close()
            p.dapiStitch = tmpStitch;
%             p.stitchDim = size(tmpStitch);
        end
        
        function p = measureChannelBackground(p, channel)
            if isempty(p.backgroundMats)
                p.backgroundMats.labels = p.stitchedScans.labels;
                p.backgroundMats.mats = cell(0, numel(p.backgroundMats.labels));
            end
            channelIdx = find(ismember(p.backgroundMats.labels, channel));
            reader = bfGetReader(p.scanFile);
            tmpStack = zeros(p.tileSize(1), p.tileSize(2), prod(p.scanDim));
            iPlane = reader.getIndex(0, channelIdx - 1, 0) + 1;
            
            for i = 1:size(tmpStack, 3)
                reader.setSeries(i-1); 
                tmpStack(:,:,i) = bfGetPlane(reader, iPlane);
            end
            backgroundMat = min(tmpStack, [], 3);
            filt = fspecial('disk',100);
            p.backgroundMats.mats{channelIdx} = uint16(imfilter(backgroundMat,filt,'replicate'));
        end
        
        function p = measureBackground(p, varargin)
            n = inputParser;
            n.addOptional('channels',  p.channels(~ismember(p.channels,{'dapi','trans'})), @(x) mustBeMember(x, p.channels))
            n.parse(varargin{:});
            
            p.backgroundMats.labels = n.Results.channels;
            p.backgroundMats.mats = cell(0, numel(p.backgroundMats.labels));
            for i = 1:numel(p.backgroundMats.labels)
                p.measureChannelBackground(p.backgroundMats.labels{i});
            end
        end
        
        function p = contrastDAPIstitch(p)
            function_scale =  @(block_struct) im2uint16(scale(block_struct.data));
            
            p.dapiStitch = blockproc(p.dapiStitch, [1000 1000], function_scale, 'BorderSize', [0 0], 'UseParallel', true);
        end
        
        function p = contrastStitchedScans(p, percentiles, scaleFactor)  %Can modify this function for different default contrast
            function_contrast =  @(block_struct) im2uint16(d2utils.percentileScaleImage(block_struct.data, percentiles, scaleFactor));
            for i = 1:numel(p.stitchedScans.stitches)
                p.stitchedScans.stitches{i} = blockproc(p.stitchedScans.stitches{i}, [1000 1000], function_contrast, 'BorderSize', [0 0], 'UseParallel', true);
            end
        end
        
        function p = contrastStitchedScanChannel(p, channel, percentiles, scaleFactor)  %Can modify this function for different default contrast
            function_contrast =  @(block_struct) im2uint16(d2utils.percentileScaleImage(block_struct.data, percentiles, scaleFactor));
            channelIdx = find(ismember(channel, p.stitchedScans.labels));
            p.stitchedScans.stitches{channelIdx} = blockproc(p.stitchedScans.stitches{channelIdx}, [1000 1000], function_contrast, 'BorderSize', [0 0], 'UseParallel', true);
        end
        
        function p = resizeStitchedScans(p)
            %Can use blockproc as below but set 'BorderSize' to [0 0]. 
%             function_resize =  @(block_struct)... 
%             imresize(block_struct.data, 1/p.resizeFactor); 
            p.smallStitchedScans.labels = p.stitchedScans.labels;
            p.smallStitchedScans.stitches = cell(0, numel(p.stitchedScans.stitches));
            for i = 1:numel(p.stitchedScans.stitches)
                p.smallStitchedScans.stitches{i} = imresize(p.stitchedScans.stitches{i}, 1/p.resizeFactor);
            end
            p.smallDapiStitch = imresize(p.dapiStitch, 1/p.resizeFactor);
%             p.smallDapiStitch = blockproc(p.dapiStitch, [5000 5000], function_resize, 'BorderSize', [0 0]);
        end
        
        function p = resizeStitchedScanChannel(p, channel)
            channelIdx = ismember(channel, p.smallStitchedScans.labels);
            p.smallStitchedScans.stitches{channelIdx} = imresize(p.stitchedScans.stitches{channelIdx}, 1/p.resizeFactor);
        end
        
        function p = resizeStitchedDapi(p)
            p.smallDapiStitch = imresize(p.dapiStitch, 1/p.resizeFactor);
        end
        
        function tmpStitch = stitchChannel(p, channel, varargin)
            n = inputParser;
            n.addRequired('channel', @(x) mustBeMember(x, p.channels));
            n.addOptional('subtractBackground', false, @islogical)
            n.parse(channel, varargin{:});
            
            channel = n.Results.channel;
            tileTable = p.tilesTable;
            tilesTmp = transpose(p.scanMatrix);
            tiles = tilesTmp(:);
            height = p.tileSize(1);
            width = p.tileSize(2);
            tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
            channelIdx = find(ismember(p.channels,channel));
            reader = bfGetReader(p.scanFile);
            iPlane = reader.getIndex(0, channelIdx - 1, 0) + 1;
            if n.Results.subtractBackground
                idx = ismember(p.backgroundMats.labels, channel);
                backgroundMat = p.backgroundMats.mats{idx};
            else
                backgroundMat = zeros(p.tileSize, 'uint16');
            end
            if logical(p.imRotation) %In case the image acquired by Elements is flipped or somehow rotated. Consider deleting. 
                for i = 1:numel(tiles) %Can this be paralelized? 
                    reader.setSeries(tiles(i)-1);
                    tmpPlane = bfGetPlane(reader, iPlane);
                    tmpPlane = tmpPlane - backgroundMat;
                    tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                    tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = rot90(tmpPlane, p.imRotation);
                end
            else
                for i = 1:numel(tiles)
                    reader.setSeries(tiles(i)-1);
                    tmpPlane = bfGetPlane(reader, iPlane);
                    tmpPlane = tmpPlane - backgroundMat;
                    tmpStitch(tileTable{tiles(i),'left'}:tileTable{tiles(i),'left'}+height-1, ...
                    tileTable{tiles(i),'top'}:tileTable{tiles(i),'top'}+width-1) = tmpPlane;
                end
            end
            reader.close()
        end

        function p = stitchChannels(p, varargin)
            n = inputParser;
            n.addOptional('channelsToStitch', p.channels(~ismember(p.channelTypes,{'dapi'})), @(x) mustBeMember(x, p.channels));
            n.parse(varargin{:});
            channelsToStitch = n.Results.channelsToStitch;
            tileTable = p.tilesTable;
            tilesTmp = transpose(p.scanMatrix);
            tiles = tilesTmp(:);
            height = p.tileSize(1);
            width = p.tileSize(2);
            stitches = cell(1,numel(channelsToStitch));
            channelIdx = find(ismember(p.channels, channelsToStitch));
            reader = bfGetReader(p.scanFile);
            for i = 1:numel(channelsToStitch) %Consider making this parfor by using bfreader memoizer
                tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
                iPlane = reader.getIndex(0, channelIdx(i) - 1, 0) + 1;
                if logical(p.imRotation) %In case the image acquired by Elements is flipped or somehow rotated. Consider deleting.
                    for ii = 1:numel(tiles)
                        reader.setSeries(tiles(ii)-1);
                        tmpPlane  = bfGetPlane(reader, iPlane);
                        tmpStitch(tileTable{tiles(ii),'left'}:tileTable{tiles(ii),'left'}+height-1, ...
                            tileTable{tiles(ii),'top'}:tileTable{tiles(ii),'top'}+width-1) = rot90(tmpPlane, p.imRotation);
                    end
                    stitches{i} = tmpStitch;
                else
                    for ii = 1:numel(tiles)
                        reader.setSeries(tiles(ii)-1);
                        tmpPlane  = bfGetPlane(reader, iPlane);
                        tmpStitch(tileTable{tiles(ii),'left'}:tileTable{tiles(ii),'left'}+height-1, ...
                            tileTable{tiles(ii),'top'}:tileTable{tiles(ii),'top'}+width-1) = tmpPlane;
                    end
                    stitches{i} = tmpStitch;
                end
            end
            reader.close()
            p.stitchedScans.labels = channelsToStitch;
            p.stitchedScans.stitches = stitches;
        end
        
        function p = stitchChannels2(p, varargin)
            %Making this a seperate methods since the background
            %subtraction step adds ~10% more time. 
            n = inputParser;
            n.addOptional('subtractBackground', false, @islogical)
            n.addOptional('channelsToStitch', p.channels(~ismember(p.channels,{'dapi','trans'})), @(x) mustBeMember(x, p.channels));
            n.parse(varargin{:});
            channelsToStitch = n.Results.channelsToStitch;
            tileTable = p.tilesTable;
            tilesTmp = transpose(p.scanMatrix);
            tiles = tilesTmp(:);
            height = p.tileSize(1);
            width = p.tileSize(2);
            stitches = cell(1,numel(channelsToStitch));
            channelIdx = find(ismember(p.channels, channelsToStitch));
            tmpBackgroundMats = zeros(height, width, numel(channelsToStitch), 'uint16');
            if n.Results.subtractBackground
                if ~isempty(p.backgroundMats)
                    for i = 1:numel(channelsToStitch)
                        idx = ismember(p.backgroundMats.labels, channelsToStitch{i});
                        tmpBackgroundMats(:,:,i) = p.backgroundMats.mats{idx};
                    end
                else
                    disp('Unable to subtract background. Please run measureBackground and try again')
                end
            end
            reader = bfGetReader(p.scanFile);
            for i = 1:numel(channelsToStitch) %Consider making this parfor by using bfreader memoizer
                tmpStitch = zeros(max(tileTable.left)+height-1,max(tileTable.top)+width-1, 'uint16');
                iPlane = reader.getIndex(0, channelIdx(i) - 1, 0) + 1;
                if logical(p.imRotation) %In case the image acquired by Elements is flipped or somehow rotated. Consider deleting.
                    for ii = 1:numel(tiles)
                        reader.setSeries(tiles(ii)-1);
                        tmpPlane  = bfGetPlane(reader, iPlane);
                        tmpPlane = tmpPlane - tmpBackgroundMats(:,:,i);
                        tmpStitch(tileTable{tiles(ii),'left'}:tileTable{tiles(ii),'left'}+height-1, ...
                            tileTable{tiles(ii),'top'}:tileTable{tiles(ii),'top'}+width-1) = rot90(tmpPlane, p.imRotation);
                    end
                    stitches{i} = tmpStitch;
                else
                    for ii = 1:numel(tiles)
                        reader.setSeries(tiles(ii)-1);
                        tmpPlane  = bfGetPlane(reader, iPlane);
                        tmpPlane = tmpPlane - tmpBackgroundMats(:,:,i);
                        tmpStitch(tileTable{tiles(ii),'left'}:tileTable{tiles(ii),'left'}+height-1, ...
                            tileTable{tiles(ii),'top'}:tileTable{tiles(ii),'top'}+width-1) = tmpPlane;
                    end
                    stitches{i} = tmpStitch;
                end
            end
            reader.close()
            p.stitchedScans.labels = channelsToStitch;
            p.stitchedScans.stitches = stitches;
        end
        
        function p = restitchChannel(p, channel, varargin)
            n = inputParser;
            n.addRequired('channel', @(x) mustBeMember(x, p.channels));
            n.addOptional('subtractBackground', false, @islogical)
            n.parse(channel, varargin{:});
            
            channel = n.Results.channel;
            channelIdx = ismember(p.stitchedScans.labels, channel);
            p.stitchedScans.stitches{channelIdx} = p.stitchChannel(channel, n.Results.subtractBackground);
        end
        
        function outIm = stitchTiles(p, rows, cols, channel, varargin) %Option to specify transformCoords - used with stitchingGUI
            if nargin == 4
                rowTransform = p.rowTransformCoords;
                colTransform = p.columnTransformCoords;
            elseif nargin > 4
                rowTransform = varargin{1};
                colTransform = varargin{2};
            end
            
            %Height and width may need to be switched for rotated images.
            %I'm not sure how tilesize changes when the image from Elements is rotated.  
            height = p.tileSize(1);
            width = p.tileSize(2);
            localScanMatrix = p.scanMatrix(rows(1):rows(2), cols(1):cols(2));
            tiles = transpose(localScanMatrix);
            tiles = tiles(:);
            topCoords = zeros(numel(tiles),1);
            leftCoords = zeros(numel(tiles),1);
            for i = 1:numel(tiles) 
                [r,c] = find(localScanMatrix == tiles(i));
                topCoords(i)  = c*colTransform(1) + r*rowTransform(1);
                leftCoords(i) = r*rowTransform(2) + c*colTransform(2);
            end
            topCoords = topCoords - min(topCoords) + 1;
            leftCoords = leftCoords - min(leftCoords) + 1;
            outIm = zeros(max(leftCoords)+height-1,max(topCoords)+width-1, 'uint16');
            reader = bfGetReader(p.scanFile);
            channelIdx = find(ismember(p.channels, channel));  
            iPlane = reader.getIndex(0, channelIdx - 1, 0) + 1;
            
            if logical(p.imRotation)
                for ii = 1:numel(tiles)
                    reader.setSeries(tiles(ii)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);
                
                    outIm(leftCoords(ii):leftCoords(ii)+height-1, ...  
                        topCoords(ii):topCoords(ii)+width-1) = rot90(tmpPlane, p.imRotation);
                end
            else
                for ii = 1:numel(tiles)
                    reader.setSeries(tiles(ii)-1);
                    tmpPlane  = bfGetPlane(reader, iPlane);
                    outIm(leftCoords(ii):leftCoords(ii)+height-1, ...
                        topCoords(ii):topCoords(ii)+width-1) = tmpPlane;
                end
            end
        end
        
        function outIm = getImageRect(p, channel, rect) %rect specified as [x y nrows ncols]
            channelIdx = ismember(p.stitchedScans.labels, channel);
            outIm = p.stitchedScans.stitches{channelIdx};
            outIm = outIm(rect(1):rect(1)+rect(3)-1, rect(2):rect(2)+rect(4)-1); %Kinda ugly, would prefer using imcrop
        end
        
        function outIm = getSmallImageRect(p, channel, rect) %rect specified as [x y nrows ncols]
            channelIdx = ismember(p.smallStitchedScans.labels, channel);
            smallRect = ceil(rect/p.resizeFactor);
            %smallRect(1:2) = max([1, 1], smallRect(1:2));
            outIm = p.smallStitchedScans.stitches{channelIdx};
            outIm = outIm(smallRect(1):smallRect(1)+smallRect(3)-1, smallRect(2):smallRect(2)+smallRect(4)-1); %Kinda ugly, would prefer using imcrop
        end
        
        function outIm = getDapiImage(p, rect) %rect specified as [x y nrows ncols]
            outIm = p.dapiStitch(rect(1):rect(1)+rect(3)-1, rect(2):rect(2)+rect(4)-1);
        end
        
        function outIm = getSmallDapiImage(p, rect) %rect specified as [x y nrows ncols]
            smallRect = ceil(rect/p.resizeFactor);
            %smallRect(1:2) = max([1, 1], smallRect(1:2));
            outIm = p.smallDapiStitch(smallRect(1):smallRect(1)+smallRect(3)-1, smallRect(2):smallRect(2)+smallRect(4)-1);
        end
        
        function [splitTiles, startPositions] = splitStitch(p, channels, varargin)
            if nargin == 2
                blockSize = min([1024 1024], p.stitchDim); %Default chunk size
            elseif narging == 3
                blockSize = varargin{1};
            end
            if all(p.scanDim == 1) %If pre-stitched scan or single tile
                rowSplit = [repmat(blockSize(1),1, floor(p.stitchDim(1)/blockSize(1))), mod(p.stitchDim(1), blockSize(1))]; 
                colSplit = [repmat(blockSize(2),1, floor(p.stitchDim(2)/blockSize(2))), mod(p.stitchDim(2), blockSize(2))];
                %Remove 0s when stitchDim is divisible by block size
                rowSplit = rowSplit(rowSplit~=0);
                colSplit = colSplit(colSplit~=0);
                splitTiles = cell(numel(rowSplit) * numel(colSplit), numel(channels));
                for i = 1:numel(channels)
                    channelIdx = ismember(p.stitchedScans.labels, channels{i});
                    tmpMat =  mat2cell(p.stitchedScans.stitches{channelIdx}, rowSplit, colSplit);
                    splitTiles(:,i) = reshape(tmpMat, 1,[]);
                end
                startPositions = combvec(1:blockSize(1):p.stitchDim(1), 1:blockSize(2):p.stitchDim(2))';
%                 startPositions = combvec(linspace(0, (nRowSplit-1)*blockSize, nRowSplit), linspace(0, (colSplit-1)*blockSize, nColSplit))';
            else %tiled scan
                splitTiles = cell(numel(p.scanMatrix), numel(channels));
                startPositions = zeros(numel(p.scanMatrix),2);
                for i = 1:numel(p.scanMatrix)
                    [r, c] = find(p.scanMatrix == i);
                    pos = p.tilesTable{i, {'left', 'top'}};
                    startPositions(i, :) = pos;
                    if c < p.scanDim(2)
                        rightTile = p.scanMatrix(r, c+1);
                        colEnd = p.tilesTable{rightTile, 'top'}-1;
                    else
                        colEnd = pos(2) + p.tileSize(2)-1;
                    end
                    if r < p.scanDim(1)
                        bottomTile = p.scanMatrix(r+1, c);
                        rowEnd = p.tilesTable{bottomTile, 'left'}-1;
                    else
                        rowEnd = pos(1) + p.tileSize(1)-1;
                    end
                    for ii = 1:numel(channels)
                        channelIdx = ismember(p.stitchedScans.labels, channels{ii});
                        splitTiles{i, ii} = p.stitchedScans.stitches{channelIdx}(pos(1):rowEnd, pos(2):colEnd);
                    end
                end
            end
        end
        
        function saveTilesTable(p)
            if ~isempty(p.tilesTable) %Not sure we need the option to specify alternative filename
                writetable(p.tilesTable, p.tilesTableName)
            else
                fprintf("tilesTable is empty. Run loadTiles and try again")
            end
        end
        
        function loadTilesTable(p)
            fprintf('Loading %s\n', p.tilesTableName);
            opts = detectImportOptions(p.tilesTableName);
            opts = setvartype(opts, 'single'); %Probably unnecessary given the typical scan size.
            p.tilesTable = readtable(p.tilesTableName, opts);
        end
        
        function loadScanSummary(p, userInputChannels)
            scanSummaryTable = d2utils.parseScanSummary(p.scanSummaryFile);
            p.scanFile = scanSummaryTable{'scanFileName',1}{:};%Could possibly check that the scanFile tilesTable exist
            p.tilesTableName = scanSummaryTable{'tilesTableName',1}{:};
            p.scanDim = str2double(split(scanSummaryTable{'scanDimensions',1})');
            p.tileSize = str2double(split(scanSummaryTable{'imageSize',1})');
            p.stitchDim = str2double(split(scanSummaryTable{'stitchDimensions',1})');
            p.snake = strcmp(scanSummaryTable{'snake',1}{:}, 'true');
            p.direction = scanSummaryTable{'scanDirection',1}{:};
            p.startPos = scanSummaryTable{'startPosition',1}{:};
            p.rowTransformCoords = str2double(split(scanSummaryTable{'rowTransform',1})');
            p.columnTransformCoords = str2double(split(scanSummaryTable{'columnTransform',1})');
            
            p.scanMatrix = d2utils.makeScanMatrix(p.scanDim, 'start', p.startPos, 'snake', p.snake,'direction', p.direction); 
            p.channels = d2utils.loadChannelsND2(p.scanFile, userInputChannels);
            
            if any(strcmp(scanSummaryTable.Properties.RowNames,'channelTypes')) %For compatibility with older version of dentist2
                p.channelTypes = split(scanSummaryTable{'channelTypes',1});
                p.checkOverChannelTypes(p.channelTypes)
            end

            %Load tiles table.
            if isfile(p.tilesTableName)
                p.loadTilesTable();
            end 
        end
        
        function saveScanSummary(p) 
            %NOTE! Avoid tabs ('\t') in your file name. Commas are OK. 
            %Should maybe check that all the necessary properties are not
            %empty. 
            
            if ~isempty(p.scanDim) && all(p.scanDim) %Using this as indicator for whether the scan is tiled or pre-stitched
                outTableArray = {p.scanFile; p.tilesTableName; num2str(p.scanDim); num2str(p.tileSize); num2str(p.stitchDim);...
                    p.startPos; p.direction; string(p.snake); strjoin(p.channels);strjoin(p.channelTypes);...
                    num2str(p.rowTransformCoords); num2str(p.columnTransformCoords)};
                outTable = cell2table(outTableArray,'RowNames', {'scanFileName', 'tilesTableName', 'scanDimensions', 'imageSize', 'stitchDimensions',...
                    'startPosition', 'scanDirection', 'snake', 'channels','channelTypes', 'rowTransform', 'columnTransform'});
            else
                outTableArray = {p.scanFile; num2str(p.stitchDim); strjoin(p.channels);strjoin(p.channelTypes)};
                outTable = cell2table(outTableArray,'RowNames', {'scanFileName', 'stitchDimensions', 'channels','channelTypes'});
            end
           
            writetable(outTable, p.scanSummaryFile, 'WriteRowNames', true, 'WriteVariableNames', false, 'QuoteStrings', false, 'Delimiter', '\t')  
        end
                
        function saveStitches(p, varargin) 
            %For the time being, saving to .mat files. 
            %At some point it may be worth modifying this method to save
            %bigTiff files for compatibility with other software. 
            %Will need to update loadStitches as well. 
            if nargin == 1 %Not sure we need the option to specify alternative filename
                outFileName = 'stitchedScans.mat';
            elseif nargin == 2
                outFileName = sprintf('%s.mat', varargin{1});
            end
            
            if ~isempty(p.dapiStitch) && ~isempty(p.stitchedScans) 
                dapi = p.dapiStitch;
                save(outFileName, 'dapi', '-v7.3')
                labels = p.stitchedScans.labels;
                save(outFileName, 'labels', '-append') 
                scans = p.stitchedScans.stitches;
                save(outFileName, 'scans', '-append')
            else
                if isempty(p.dapiStitch)
                    fprintf("dapiStitch is empty. Run stitchDAPI and try again\n")
                end
                if isempty(p.stitchedScans)
                    fprintf("stitchedScans is empty. Run stitchChannels and try again\n")
                end
            end
       end
       
       function p = loadStitches(p)
           tmpMat = matfile('stitchedScans.mat');
           p.dapiStitch = tmpMat.dapi;
           p.stitchedScans = struct;
           p.stitchedScans.labels = tmpMat.labels;
           p.stitchedScans.stitches = tmpMat.scans;
       end
       
       function p = reloadChannelStitch(p, channel)
           tmpMat = matfile('stitchedScans.mat');
           labels = tmpMat.labels;
           idx = find(ismember(channel, labels));
           tmpStitch = tmpMat.scans(1,idx);
           p.stitchedScans.stitches{idx} = tmpStitch{:};
       end
        
    end
         
end