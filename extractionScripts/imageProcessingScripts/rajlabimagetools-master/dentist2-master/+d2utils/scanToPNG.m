function [] = scanToPNG(varargin)
    n = inputParser;
    n.addParameter('scanSummary', 'scanSummary.txt', @ischar);
    n.addParameter('preStitchedScan', '', @ischar);
    n.addParameter('stitchFile', 'stitchedScans.mat', @ischar);
    n.addParameter('channel', '', @ischar);
    n.addParameter('outDir', 'cp', @ischar);
    n.addParameter('outPrefix', '', @ischar);
%     n.addParameter('splitFactor', [1 1], @(x)validateattributes(x,{'numeric'}, {'size', [1,2], '>', 0}));
    n.addParameter('contrastPercentile', [1 99], @(x)validateattributes(x,{'numeric'}, {'size', [1,2], 'nonnegative'}));
    n.addParameter('contrastScale', [1 1], @(x)validateattributes(x,{'numeric'}, {'size', [1,2], 'nonnegative'}));
    n.addParameter('resizeFactor', 1, @(x)validateattributes(x,{'numeric'}, {'scalar', '>', 0}));

    n.parse(varargin{:});
    if isempty(n.Results.preStitchedScan)
        if isfile(n.Results.scanSummary)
            scanObj = scanObject('scanSummary', n.Results.scanSummary);
        else
            fprintf('Unable to detect %s in your current directory.\n. Make sure to run the d2StitchingGUI before running this function.\n You may also want to check your path and %s and try again. ', n.Results.scanSummary, n.Results.scanSummary)
            return
        end

        %Check if stitches have been saved. If not, use scanObj methods to stitch.
        if isempty(scanObj.tilesTable)
            disp('The scan object does not contain a tiles table. Creating a new tiles table.')
            scanObj.loadTiles();
            scanObj.savetilesTable();
        end

        if isfile('stitchedScans.mat')
            fprintf('Loading stitched scan.\n')
            %Note that this doesn't check if channel exists in the stitchedScans.mat file.
            tmpStitch = d2utils.loadStitchD2(n.Results.stitchFile, n.Results.channel);
        else
            fprintf('Stitching %s channel.', n.Results.channel)
            tmpStitch = scanObj.stitchChannel(n.Results.channel);
        end
    else
        channels = d2utils.readND2Channels(n.Results.preStitchedScan);
        if ismember(n.Results.channel, channels)
            fprintf('Loading %s image from %s.\n', n.Results.channel, n.Results.preStitchedScan)
            channelIdx = find(ismember(channels, n.Results.channel));
            tmpStitch = d2utils.ND2tileReader(n.Results.preStitchedScan, 'channel', channelIdx);
        else
            fprintf('Unable to find channel "%s" in scan %s.\nPlease specify one of these available channels: %s.\n', n.Results.channel, n.Results.preStitchedScan, string(join(channels, ', ')))
            return
        end
        
    end
    %Consider adding section to split scan 
    %Resize stitch
    if ~(n.Results.resizeFactor == 1)
        tmpStitch = imresize(tmpStitch, 1/n.Results.resizeFactor);
    end
    %Contrast stitch
    if all(n.Results.contrastPercentile > 0) && all(n.Results.contrastPercentile > 0)
        tmpStitch = im2uint16(d2utils.percentileScaleImage(tmpStitch, n.Results.contrastPercentile, n.Results.contrastScale));
    end
    %Save stitch
    if isempty(n.Results.outPrefix)
        outPref = n.Results.channel;
    else
        outPref = n.Results.outPrefix;
    end
    
    if ~exist(n.Results.outDir, 'dir')
        mkdir(n.Results.outDir)
    end
    imwrite(tmpStitch, fullfile(n.Results.outDir, sprintf('%s.png', outPref)))
end