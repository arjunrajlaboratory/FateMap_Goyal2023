function mask = makeDAPImask2(dapiIm, varargin)
    
    p = inputParser;
    p.addRequired('dapiIm', @(x)validateattributes(x,{'numeric'}, {'2d'}));
    p.addParameter('sensitivity', 0.1, @(x)validateattributes(x,{'numeric'}, {'scalar', '>=',0,'<=',1.0})); 
    p.addParameter('strelRad', 0, @(x)validateattributes(x,{'double'}, {'scalar', '>=',0})); 
    
    p.parse(dapiIm, varargin{:});
    
    if p.Results.strelRad > 0
        se = strel('disk', p.Results.strelRad);
        tmpOpen = imopen(p.Results.dapiIm,se);
        tmpDilate = imdilate(p.Results.dapiIm,se);
        T = adaptthresh(tmpDilate, p.Results.sensitivity,'ForegroundPolarity','bright');
        mask = imbinarize(scale(tmpOpen),T);
    else
        T = adaptthresh(p.Results.dapiIm, p.Results.sensitivity,'ForegroundPolarity','bright');
        mask = imbinarize(p.Results.dapiIm,T);
    end
    
end

