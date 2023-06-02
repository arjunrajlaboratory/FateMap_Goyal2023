function outColors = expressionToColors(expressionVector, palette)
    switch(palette)
        case 'BuYlRd'
            sampleColors = [0.1922 0.2118 0.5843; 0.2706 0.4588 0.7059; 0.45 0.68 0.82;...
                0.67 0.85 0.91; 0.88 0.95 0.97; 0.99 0.88 0.56; 0.99 0.68 0.38;...
                0.96 0.43 0.26; 0.84 0.19 0.15; 0.65 0 0.15];
        case 'YlOrRd'
            sampleColors = [1 0.85 0.46; 1 0.7 0.3; 0.99 0.55 0.23; 0.99 0.31 0.16; 0.89 0.1 0.11; 0.7 0 0.15];
        case 'GrBu'
            sampleColors = [0.8275 0.8275 0.8275; 0 0 1];
        case 'BuGnYlRd'
            sampleColors = [0 0 1; 0 0.75 1; 0.165 1 0.189; 1 0.85 0; 1 0.66 0.04; 1 0 0];
        otherwise
            sampleColors = [0.1922 0.2118 0.5843; 0.2706 0.4588 0.7059; 0.45 0.68 0.82;...
                0.67 0.85 0.91; 0.88 0.95 0.97; 0.99 0.88 0.56; 0.99 0.68 0.38;...
                0.96 0.43 0.26; 0.84 0.19 0.15; 0.65 0 0.15];
    end
    
    maxCount = max(expressionVector);
    if maxCount == 0
        outColors = repmat(sampleColors(1), numel(expressionVector), 1);
    elseif maxCount < height(sampleColors)
        outColors = interp1(0:maxCount, sampleColors(1:maxCount+1), expressionVector);
    else
        outColors = interp1(round(linspace(0, maxCount, height(sampleColors))), sampleColors, expressionVector);
    end
    outColors = single(outColors);
end