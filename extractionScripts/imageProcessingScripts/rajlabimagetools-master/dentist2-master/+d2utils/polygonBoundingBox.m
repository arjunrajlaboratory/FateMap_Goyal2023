function outRect = polygonBoundingBox(polygon)
    mins = min(polygon);
    maxs = max(polygon);
    outRect = [floor(mins), ceil(maxs-mins)];
end