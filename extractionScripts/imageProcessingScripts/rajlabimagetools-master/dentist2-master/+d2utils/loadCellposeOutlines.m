function polyArray = loadCellposeOutlines(inFileName, varargin)
    p = inputParser;
    p.addRequired('inFileName', @ischar);
    p.addParameter('position', [1 1],@(x)validateattributes(x,{'numeric'}, {'size', [1, 2], '>', 0}));
    p.addParameter('scaleFactor', 1, @(x)validateattributes(x,{'numeric'}, {'scalar', '>', 0}));
    p.parse(inFileName, varargin{:});
    
    inFileName = p.Results.inFileName;
    scaleFactor = p.Results.scaleFactor;
    position = ceil(p.Results.position/scaleFactor);
    
    tmpFile = dir(inFileName);
    if logical(tmpFile.bytes) %trouble with empty files. Should try improving how files is read.
        f = fopen(inFileName);
        tmpTable = cell2mat(textscan(f, '', 'Delimiter', ','));
        fclose(f);
        
        tmpTable = array2table(tmpTable);
        polyArray = rowfun(@(x) reshape(x, 2, [])', tmpTable, 'SeparateInputs', false, 'OutputFormat','cell');
        polyArray = cellfun(@(x) x(~isnan(x(:,1)), :) + position,polyArray, 'UniformOutput', false);
    else
        polyArray = {};
    end
end
