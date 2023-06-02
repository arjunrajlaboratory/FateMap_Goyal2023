function outRect = expandView2x(inRect, maxDim)
    startPos = max([1, 1], [inRect(1) - (inRect(3)/2), inRect(2) - (inRect(4)/2)]); %Avoid rect outside range of scan. 
    startPos = min(maxDim- inRect(3:4) +1,round(startPos));
    scanSize = min(maxDim - startPos+1, 2 * inRect(3:4));
    outRect = [startPos, scanSize];
end