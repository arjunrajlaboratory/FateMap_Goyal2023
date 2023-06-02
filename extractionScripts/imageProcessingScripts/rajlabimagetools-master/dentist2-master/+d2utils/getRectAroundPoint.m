function outRect = getRectAroundPoint(point, height, width, maxDimensions)
    startPos = max([1, 1], point - [height/2 width/2] + 1); %Avoid rect outside range of scan. 
    startPos = min(maxDimensions - [height width],startPos);
    outRect = single(round([startPos, height, width]));
end