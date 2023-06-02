function outMask = polyshape2mask(polyIn, shift, sz, varargin)
    p = inputParser;
    p.addRequired('polyIn', @(x)validateattributes(x,{'polyshape'}, {'scalar'}));
    p.addRequired('shift', @(x)validateattributes(x,{'numeric'}, {'size', [1,2]}));
    p.addRequired('sz', @(x)validateattributes(x,{'numeric'}, {'size', [1,2]}));
    p.addParameter('flip', false, @islogical);
    p.parse(polyIn, shift, sz, varargin{:});
    
    polyHoles = holes(p.Results.polyIn); %Not necessary if we don't allows for holes... 
    polyIn = rmholes(p.Results.polyIn);
    if p.Results.flip
        tmpVertices = fliplr(polyIn.Vertices - p.Results.shift);
    else
        tmpVertices = polyIn.Vertices - p.Results.shift;
    end
    sz = p.Results.sz;
    outMask = poly2mask(tmpVertices(:,1), tmpVertices(:,2), sz(1),  sz(2));
    if ~isempty(polyHoles)
        if p.Results.flip
            for i = 1:numel(polyHoles)
                tmpHoleVertices = fliplr(polyHoles(i).Vertices - p.Results.shift);
                tmpHoleMask = poly2mask(tmpHoleVertices(:,1), tmpHoleVertices(:,2), sz(1),  sz(2));
                outMask = outMask & ~tmpHoleMask;
            end
        else
            for i = 1:numel(polyHoles)
                tmpHoleVertices = polyHoles(i).Vertices - p.Results.shift;
                tmpHoleMask = poly2mask(tmpHoleVertices(:,1), tmpHoleVertices(:,2), sz(1),  sz(2));
                outMask = outMask & ~tmpHoleMask;
            end
        end
    end
end