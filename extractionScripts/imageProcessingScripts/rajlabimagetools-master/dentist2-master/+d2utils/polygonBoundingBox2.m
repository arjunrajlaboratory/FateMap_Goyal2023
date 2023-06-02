function outRect = polygonBoundingBox2(polygon)
    mins = min(polygon);
    maxs = max(polygon);
    outRect = [mins-1, maxs-mins+1];
end