function outRect = polyshapeBoundingBox(polyin, varargin)
    p = inputParser;
    p.addRequired('polyin', @(x)validateattributes(x,{'polyshape'}, {'vector'}));
    p.addOptional('maxDim', [inf, inf], @(x) validateattributes(x,{'numeric'},  {'size', [1,2]}));
    p.parse(polyin, varargin{:});
    
    [xlim, ylim] = boundingbox(p.Results.polyin);
    if isempty(xlim) %Assuming that if xlim is not empty, then ylim is not empty
        outRect = zeros(1,4);
    else
        start = floor(max([xlim(1), ylim(1)], [1, 1]));
        sz = ceil([diff(xlim), diff(ylim)]);
        sz = min(sz, p.Results.maxDim-start);
        outRect = [start, sz];
    end
    
end