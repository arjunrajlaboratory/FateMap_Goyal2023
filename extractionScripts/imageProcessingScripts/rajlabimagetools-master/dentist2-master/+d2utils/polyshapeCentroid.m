function out = polyshapeCentroid(polyIn)
    [x, y] = centroid(polyIn, 1:polyIn.NumRegions);
    out = [x, y];
end