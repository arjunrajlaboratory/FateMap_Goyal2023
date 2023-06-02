function  coordsOut = poly2centroid(polyIn)
    [xOut, yOut] = centroid(polyshape(polyIn));
    coordsOut = [xOut, yOut];
end