function [localCoordsX, localCoordsY] = globalToLocalCoords(rect,globalCoordsX,globalCoordsY)
    localCoordsX = globalCoordsX - rect(1) + 1;
    localCoordsY = globalCoordsY - rect(2) + 1;

end