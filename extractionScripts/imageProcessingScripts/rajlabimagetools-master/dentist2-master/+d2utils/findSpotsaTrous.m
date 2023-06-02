function [X, Y, spotInt] = findSpotsaTrous(img, varargin)
    p = inputParser;
    p.addRequired('img', @(x)validateattributes(x,{'numeric'}, {'2d'}));
    p.addParameter('numLevels', 3, @(x)validateattributes(x,{'numeric'}, {'scaler', '>',0})); 
    p.addParameter('sigma', 2, @(x)validateattributes(x,{'numeric'}, {'scaler', '>',0}));
    p.addParameter('shift', [0, 0], @(x)validateattributes(x,{'numeric'}, {'size', [1,2]}))
    p.addParameter('threshFactor', 2, @(x)validateattributes(x,{'numeric'}, {'scalar', '>', 0}))

    p.parse(img, varargin{:});
    
    img = p.Results.img;
    numLevels = p.Results.numLevels;
    sigma = p.Results.sigma;
    shift = p.Results.shift;
    
    [aTrous, ~] = aTrousWaveletTransform(img,'numLevels',numLevels,'sigma',sigma);
    imgAT = sum(aTrous,3);
    bw = imregionalmax(imgAT);
    regionalMaxValues = imgAT(bw);
    regionalMaxIndices = find(bw);

    [regionalMaxValues,I] = sort(regionalMaxValues,'ascend');

    regionalMaxIndices = regionalMaxIndices(I);
    %Auto threshold
    [threshold] = imregmaxThresh(regionalMaxValues);
    if isempty(threshold)
        threshold = max(regionalMaxValues) + 1; %beyond max
    end
    
    minRegionalMax = floor(threshold/p.Results.threshFactor);
    spotInds = regionalMaxIndices(regionalMaxValues>minRegionalMax);
    spotInt = regionalMaxValues(regionalMaxValues>minRegionalMax);

    [X, Y] = ind2sub(size(bw),spotInds);  % convert 1D to 2D
    if any(shift)
        X = X + shift(1);
        Y = Y + shift(2);
    end
end