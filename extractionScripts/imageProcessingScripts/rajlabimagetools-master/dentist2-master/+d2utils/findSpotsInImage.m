function [X, Y, intensities] = findSpotsInImage(img, percentileToKeep, varargin)
    p = inputParser;
    p.addRequired('img', @(x)validateattributes(x,{'numeric'}, {'2d'}));
    p.addRequired('percentileToKeep', @(x)validateattributes(x,{'numeric'}, {'scalar', '>=',1,'<=',100})); 
    p.addParameter('filterSize', 20, @(x)validateattributes(x,{'numeric'}, {'scalar','integer', '>',0})); 
    p.addParameter('sigma', 2, @(x)validateattributes(x,{'numeric'}, {'scaler', '>',0})); 
    p.addParameter('filter', [], @(x)validateattributes(x,{'double'}, {'2d'})); 
    p.addParameter('shift', [0, 0], @(x)validateattributes(x,{'numeric'}, {'size', [1,2]}))
    p.addParameter('connectivity', 8, @(x)validateattributes(x,{'numeric'}))
    p.parse(img, percentileToKeep, varargin{:});
    
    img = p.Results.img;
    percentileToKeep = p.Results.percentileToKeep;
    filterSize = p.Results.filterSize;
    sigma = p.Results.sigma;
    shift = p.Results.shift;
    conn = p.Results.connectivity;
    
    if isempty(p.Results.filter)
        theFilter = -fspecial('log',filterSize,sigma);
    else
        theFilter = p.Results.filter;
    end
    
    filt = imfilter(im2single(img),theFilter,'replicate');
    irm = imregionalmax(filt, conn);
    tempSpots = filt(irm)';
    thresh = prctile(tempSpots,percentileToKeep); %I wonder if we should put this on log scale
    filt = filt.*single(irm);
%     filt(filt < thresh) = 0;
%     filt = im2uint16(filt);
    idx = filt >= thresh;
    intensities = im2uint16(filt(idx));
    [X,Y] = ind2sub(size(filt),find(idx));
    if any(shift)
        X = X + shift(1);
        Y = Y + shift(2);
    end
end