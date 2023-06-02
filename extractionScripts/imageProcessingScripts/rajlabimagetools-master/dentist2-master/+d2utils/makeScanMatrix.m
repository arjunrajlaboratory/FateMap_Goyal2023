    function scanMatrix = makeScanMatrix(dimensions, varargin)

    p = inputParser;

    p.addRequired('dimensions', @(x)validateattributes(x,{'numeric'}, {'size',[1 2]})); %[rows columns]

    p.addParameter('start', 'top left', @(x) assert(ismember(lower(x), {'top left', 'top right', ...
        'bottom left', 'bottom right'}), 'Options are: "top left", "top right", "bottom left", "bottom right"'));
    p.addParameter('snake', true, @islogical);
    p.addParameter('direction', 'horizontal', @(x) assert(ismember(lower(x), {'horizontal', 'vertical'}),...
        'Options are: "horizontal" or "vertical"'));
    p.addParameter('outFile', '', @(x) assert((ischar(x) & endsWith(x, '.csv')),...
        'Specify .csv filename to save scanMatrix. e.g "scanMatrix.csv"'));

    p.parse(dimensions, varargin{:});

    scanDim = p.Results.dimensions;
    start = lower(p.Results.start);
    direction = lower(p.Results.direction);
    
    %scanMatrix = vec2mat(1:scanDim(1)*scanDim(2), scanDim(2)); % This requires communications toolbox
    scanMatrix=reshape(1:scanDim(1)*scanDim(2), scanDim(2),scanDim(1))'; % doesn't require communications toolbox
    
    if strcmp(direction, 'horizontal')

        if p.Results.snake
            if ismember(start, {'top right', 'bottom right'})
                for i = 1:2:scanDim(1)
                    scanMatrix(i, :) = fliplr(scanMatrix(i, :));
                end
          elseif ismember(start, {'top left', 'bottom left'})
              for i = 2:2:scanDim(1)
                scanMatrix(i, :) = fliplr(scanMatrix(i, :));
              end  
            end
        else
            if ismember(start, {'top right', 'bottom right'})
                scanMatrix = fliplr(scanMatrix);
            end
        end

        if ismember(start, {'bottom left', 'bottom right'})
            scanMatrix = flipud(scanMatrix);
        end

    elseif strcmp(direction, 'vertical')
        scanMatrix = scanMatrix.';

        if p.Results.snake
            if ismember(start, {'bottom right', 'bottom left'})
                for i = 1:2:scanDim(2)
                    scanMatrix(:, i) = flipud(scanMatrix(:, i));
                end
          elseif ismember(start, {'top right', 'top left'})
              for i = 2:2:scanDim(2)
                scanMatrix(:, i) = flipud(scanMatrix(:, i));
              end  
            end
        else
            if ismember(start, {'bottom right', 'bottom left'})
                scanMatrix = flipud(scanMatrix);
            end
        end

        if ismember(start, {'top right', 'bottom right'})
            scanMatrix = fliplr(scanMatrix);
        end

    end
    
    if ~isempty(p.Results.outFile)
        writematrix(scanMatrix, p.Results.outFile)
    end
    
end
