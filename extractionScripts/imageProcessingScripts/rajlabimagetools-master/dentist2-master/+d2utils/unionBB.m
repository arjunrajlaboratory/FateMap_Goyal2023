function outBB = unionBB(matBB)
%matBB should be a nx4 matrix where each row represents a bounding box:
%[x, y, height, width];

tmpStart = min(matBB(:,1:2), [], 1);
tmpEnd = matBB(:,1:2) + matBB(:,3:4);
tmpEnd = max(tmpEnd, [], 1);
sz = tmpEnd - tmpStart;
outBB = [tmpStart, sz];

end