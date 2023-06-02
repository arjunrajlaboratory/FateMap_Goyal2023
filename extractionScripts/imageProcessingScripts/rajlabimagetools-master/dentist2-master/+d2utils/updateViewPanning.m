function outRect = updateViewPanning(inRect, displacement, maxDim)
    startPos = max([1, 1], [inRect(1) + displacement(2), inRect(2) + displacement(1)]); %Avoid rect outside range of scan. 
    startPos = min(maxDim- inRect(3:4),round(startPos));
    outRect = single([startPos, inRect(3:4)]);
end