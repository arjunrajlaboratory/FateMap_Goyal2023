function polygonVector = parseCellposeOutlines(inFileName, varargin)
    p = inputParser;
    p.addRequired('inFileName', @ischar);
    p.addParameter('position', [1 1],@(x)validateattributes(x,{'numeric'}, {'size', [1, 2], '>', 0}));
    p.addParameter('scaleFactor', 1, @(x)validateattributes(x,{'numeric'}, {'vector', '>', 0}));
    p.addParameter('flip', false, @islogical);
    p.parse(inFileName, varargin{:});
    
    inFileName = p.Results.inFileName;
    scaleFactor = p.Results.scaleFactor;
    position = ceil(p.Results.position./scaleFactor);
    
    tmpFile = dir(inFileName);
    if logical(tmpFile.bytes) %trouble with empty files. Should try improving how files is read.
        f = fopen(inFileName);
        tmpTable = cell2mat(textscan(f, '', 'Delimiter', ','));
        fclose(f);
        
        tmpTable = array2table(tmpTable);
        tmpTable = rowfun(@(x) reshape(x, 2, [])', tmpTable, 'SeparateInputs', false, 'OutputFormat','cell');
        if p.Results.flip
            tmpTable = cellfun(@(x) fliplr(x),tmpTable, 'UniformOutput', false);
        end
%         tmpTable = cellfun(@(x) x + position,tmpTable, 'UniformOutput', false);
        tmpTable = cellfun(@(x) x(~isnan(x(:,1)), :) + position,tmpTable, 'UniformOutput', false);
        tmpTable = cellfun(@(x) polyshape(x),tmpTable, 'UniformOutput', false);
        polygonVector = [tmpTable{:}];
    else
        polygonVector = [];
    end
end

