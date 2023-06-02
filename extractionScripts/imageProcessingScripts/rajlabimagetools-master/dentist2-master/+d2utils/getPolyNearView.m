function [polyArrayOut, polyNearIdx] = getPolyNearView(polyArrayIn, startPos, size)
    startPos = startPos-size/4;
    startPos = max([1, 1], startPos);
    boundingRect = polyshape([startPos(1), startPos(1)+ size(1), startPos(1)+ size(1),  startPos(1)],...
        [startPos(2), startPos(2), startPos(2)+ size(2), startPos(2)+ size(2)]);
    polyNearIdx = overlaps(boundingRect, polyArrayIn);
    polyArrayOut = polyArrayIn(polyNearIdx);
end
