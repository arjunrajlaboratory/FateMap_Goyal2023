function polyOut = largestPolyRegion(polyIn)
    tmpPoly = regions(sortregions(polyIn, 'area','descend'));
    polyOut = tmpPoly(1);
end