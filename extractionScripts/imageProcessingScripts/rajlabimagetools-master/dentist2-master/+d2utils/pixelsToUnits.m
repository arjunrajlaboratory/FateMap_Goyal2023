function units = pixelsToUnits(pixelMat, figureSize)
    %Pixel mat is [left bottom width height] and figureSize is [width
    %height]
    units(1) = pixelMat(1)/figureSize(1);
    units(2) = pixelMat(2)/figureSize(2);
    units(3) = pixelMat(3)/figureSize(1);
    units(4) = pixelMat(4)/figureSize(2);
end