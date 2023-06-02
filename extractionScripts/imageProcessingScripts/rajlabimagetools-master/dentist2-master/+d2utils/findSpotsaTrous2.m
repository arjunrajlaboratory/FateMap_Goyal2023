function [X, Y, spotInt] = findSpotsaTrous2(img,varargin)
    p = inputParser;
    p.addRequired('img', @(x)validateattributes(x,{'numeric'}, {'2d'}));
    p.addParameter('sigma',2,@(x)validateattributes(x,{'numeric'}, {'scalar','>',0}))
    p.addParameter('numLevels', 3, @(x)validateattributes(x,{'numeric'}, {'scalar', '>',0}))
    p.addParameter('shift', [0, 0], @(x)validateattributes(x,{'numeric'}, {'size', [1,2]}))
    p.addParameter('minRegionalMax',[],@(x) isempty(x) || all([isnumeric(x), isscalar(x), x>=0]))
    p.addParameter('threshFactor', 2, @(x)validateattributes(x,{'numeric'}, {'scalar', '>', 0}))

    p.parse(img, varargin{:});
    
    img = p.Results.img;
    sigma = p.Results.sigma;
    numLevels = p.Results.numLevels;
    shift = p.Results.shift;
    minRegionalMax=p.Results.minRegionalMax;
    threshFactor=p.Results.threshFactor;
    
    
    [aTrous, ~] = aTrousWaveletTransform(img,'numLevels',numLevels,'sigma',sigma);
    imgAT = sum(aTrous,3);
    bw = imregionalmax(imgAT);
    regionalMaxValues = imgAT(bw);
    regionalMaxIndices = find(bw);

    [regionalMaxValues,I] = sort(regionalMaxValues,'ascend');

    regionalMaxIndices = regionalMaxIndices(I);
    
    if isempty(minRegionalMax)
        %Auto threshold
        [autoThreshold] = imregmaxThresh(regionalMaxValues);
        
        if isempty(autoThreshold)
            autoThreshold = max(regionalMaxValues) + 1; %beyond max
        end
        minRegionalMax = floor(autoThreshold/threshFactor);
    end
        
    spotInds = regionalMaxIndices(regionalMaxValues>minRegionalMax);
    spotInt = regionalMaxValues(regionalMaxValues>minRegionalMax);

    [X, Y] = ind2sub(size(bw),spotInds);  % convert 1D to 2D
    if any(shift)
        X = X + shift(1);
        Y = Y + shift(2);
    end
end