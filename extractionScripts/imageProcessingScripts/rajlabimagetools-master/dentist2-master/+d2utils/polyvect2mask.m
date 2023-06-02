function outMask = polyvect2mask(sz, varargin)
    p = inputParser;
    p.addRequired('sz', @(x)validateattributes(x,{'numeric'}, {'size', [1,2]}));
    p.addOptional('polyVect', polyshape(), @(x) validateattributes(x,{'polyshape'}, {'vector'}));
    p.addOptional('shift', [0 0], @(x)validateattributes(x,{'numeric'}, {'size', [1,2]}));
    p.addParameter('flip', false, @islogical);
    p.parse(sz, varargin{:});
    
    sz = p.Results.sz;
    polyVect = p.Results.polyVect;
    outMask = false(sz);
    for i = 1:numel(polyVect)
        tmpMask = d2utils.polyshape2mask(polyVect(i), p.Results.shift, sz,  'flip', p.Results.flip);
        outMask = or(outMask, tmpMask);
    end
end