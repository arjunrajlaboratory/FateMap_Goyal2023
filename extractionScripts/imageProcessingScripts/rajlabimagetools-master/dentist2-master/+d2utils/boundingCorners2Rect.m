function outRect = boundingCorners2Rect(BB)
    mins = min(BB, [], 1);
    maxs = max(BB, [], 1);
    outRect = [mins, maxs-mins];
end