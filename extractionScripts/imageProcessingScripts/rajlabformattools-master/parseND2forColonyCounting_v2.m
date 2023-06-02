function [] = parseND2forColonyCounting_v2(inFile, dimensions, varargin)
% Script to take scan acquired using Nikon elements and convert to
% tiff format compatible with Rajlab image tools. Require bfmatlab code is in
% your matlab path for reading nd2 files. The bfmatlab code can be
% downloaded here https://docs.openmicroscopy.org/bio-formats/5.7.3/developers/matlab-dev.html
% Includes options to specify scanning pattern and whether or not to rotate
% each individual image tile (occasionally a problem with ELements). 

% Examples of usage:
% stitchNd2ForColonyCounting('20201115_184113_906__WellA2_ChannelDAPI_Seq0003.nd2', [29, 29], 'outDir', 'WellA2')
% stitchNd2ForColonyCounting('20201115_184113_906__WellA2_ChannelDAPI_Seq0003.nd2', [29, 29], 'outDir', 'WellA2','start', 'top left', 'snake', true, 'direction', 'horizontal')

p = inputParser;

p.addRequired('inFile', @ischar);
p.addRequired('dimensions', @(x)validateattributes(x,{'numeric'}, {'size',[1 2]}));

p.addParameter('outDir', '', @ischar);
p.addParameter('start', 'top left', @(x) assert(ismember(lower(x), {'top left', 'top right', ...
    'bottom left', 'bottom right'}), 'Options are: "top left", "top right", "bottom left", "bottom right"'));
p.addParameter('snake', true, @islogical);
p.addParameter('direction', 'horizontal', @(x) assert(ismember(lower(x), {'horizontal', 'vertical'}),...
    'Options are: "horizontal" or "vertical"'));
p.addParameter('rotate', 0, @isnumeric);

p.parse(inFile, dimensions, varargin{:});

scanDim = p.Results.dimensions;

if ~isempty(p.Results.outDir)
    outDir = p.Results.outDir;
else
    outDir = 'splitScan';
end

if ~isdir(outDir)
    mkdir(outDir)
end

%Create matrix of how the scan was acquired. This can be updated as a
scanMatrix = vec2mat(1:scanDim(1)*scanDim(2), scanDim(2));

if strcmp(p.Results.direction, 'horizontal')
    
    if p.Results.snake
        if ismember(p.Results.start, {'top right', 'bottom right'})
            for i = 1:2:scanDim(1)
                scanMatrix(i, :) = fliplr(scanMatrix(i, :));
            end
      elseif ismember(p.Results.start, {'top left', 'bottom left'})
          for i = 2:2:scanDim(1)
            scanMatrix(i, :) = fliplr(scanMatrix(i, :));
          end  
        end
    else
        if ismember(p.Results.start, {'top right', 'bottom right'})
            scanMatrix = fliplr(scanMatrix);
        end
    end

    if ismember(p.Results.start, {'bottom left', 'bottom right'})
        scanMatrix = flipud(scanMatrix);
    end
    
elseif strcmp(p.Results.direction, 'vertical')
    scanMatrix = scanMatrix.';
    
    if p.Results.snake
        if ismember(p.Results.start, {'bottom right', 'bottom left'})
            for i = 1:2:scanDim(2)
                scanMatrix(:, i) = flipud(scanMatrix(:, i));
            end
      elseif ismember(p.Results.start, {'top right', 'top left'})
          for i = 2:2:scanDim(2)
            scanMatrix(:, i) = flipud(scanMatrix(:, i));
          end  
        end
    else
        if ismember(p.Results.start, {'bottom right', 'bottom left'})
            scanMatrix = flipud(scanMatrix);
        end
    end
    
    if ismember(p.Results.start, {'top right', 'bottom right'})
        scanMatrix = fliplr(scanMatrix);
    end
    
end

tiles = scanMatrix(:);

% Read nd2 file and write to tiff in order expected by colonycounting_v2
reader = bfGetReader(p.Results.inFile); 

omeMeta = reader.getMetadataStore();
wavelengths = omeMeta.getPixelsSizeC(0).getValue(); % number of wavelength channels

for ii = 1:numel(tiles)

    for iii = 1:wavelengths

    reader.setSeries(tiles(ii)-1);
    iPlane = reader.getIndex(0, iii - 1, 0) + 1;
    tmpPlane  = bfGetPlane(reader, iPlane);

    if ~p.Results.rotate == 0
       tmpPlane = imrotate(tmpPlane, p.Results.rotate);
    end

    imwrite(tmpPlane, fullfile(outDir, sprintf('Scan001_w%d_s%d_t1.TIF', iii, ii)))

    end 

end
    


end

