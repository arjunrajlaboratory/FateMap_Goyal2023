function scanMatrix = makeScanMatrix(scanObject)
    scanDim = scanObject.dimensions;
    orientation = scanObject.orientation;

    %Create matrix of how the scan was acquired. 
    scanMatrix = vec2mat(1:scanDim(1)*scanDim(2), scanDim(2));

    if strcmp(orientation.direction, 'horizontal')

        if orientation.snake
            if ismember(orientation.start, {'top right', 'bottom right'})
                for i = 1:2:scanDim(1)
                    scanMatrix(i, :) = fliplr(scanMatrix(i, :));
                end
          elseif ismember(orientation.start, {'top left', 'bottom left'})
              for i = 2:2:scanDim(1)
                scanMatrix(i, :) = fliplr(scanMatrix(i, :));
              end  
            end
        else
            if ismember(orientation.start, {'top right', 'bottom right'})
                scanMatrix = fliplr(scanMatrix);
            end
        end

        if ismember(orientation.start, {'bottom left', 'bottom right'})
            scanMatrix = flipud(scanMatrix);
        end

    elseif strcmp(orientation.direction, 'vertical')
        scanMatrix = scanMatrix.';

        if orientation.snake
            if ismember(orientation.start, {'bottom right', 'bottom left'})
                for i = 1:2:scanDim(2)
                    scanMatrix(:, i) = flipud(scanMatrix(:, i));
                end
          elseif ismember(orientation.start, {'top right', 'top left'})
              for i = 2:2:scanDim(2)
                scanMatrix(:, i) = flipud(scanMatrix(:, i));
              end  
            end
        else
            if ismember(orientation.start, {'bottom right', 'bottom left'})
                scanMatrix = flipud(scanMatrix);
            end
        end

        if ismember(orientation.start, {'top right', 'bottom right'})
            scanMatrix = fliplr(scanMatrix);
        end

    end

end