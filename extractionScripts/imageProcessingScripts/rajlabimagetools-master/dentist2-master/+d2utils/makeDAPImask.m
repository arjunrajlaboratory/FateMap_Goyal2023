function mask = makeDAPImask(dapiIm, varargin)
    
    p = inputParser;
    p.addRequired('dapiIm', @(x)validateattributes(x,{'numeric'}, {'2d'}));
    p.addParameter('sensitivity', 0.1, @(x)validateattributes(x,{'numeric'}, {'scalar', '>=',0,'<=',1.0}));    
    
    p.parse(dapiIm, varargin{:});
    
    T = adaptthresh(p.Results.dapiIm, p.Results.sensitivity,'ForegroundPolarity','bright');
    mask = imbinarize(p.Results.dapiIm,T);
end

