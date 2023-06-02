function [globalCoordsX,globalCoordsY] = localToGlobalCoords(rect,localCoordsX,localCoordsY)
    globalCoordsX = localCoordsX + rect(1) - 1;
    globalCoordsY = localCoordsY + rect(2) - 1;
end