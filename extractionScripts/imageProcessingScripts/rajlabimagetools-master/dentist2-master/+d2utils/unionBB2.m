function outBB = unionBB2(matBB, varargin)
%matBB should be a nx4 matrix where each row represents a bounding box:
%[x, y, height, width];
%varargin{1} is boundary added to [x, y] 
if numel(varargin) > 0
    boundary = varargin{1};
else
    boundary = [0 0];
end
tmpStart = min(matBB(:,1:2), [], 1) - boundary;
tmpEnd = matBB(:,1:2) + matBB(:,3:4) + boundary;
tmpEnd = max(tmpEnd, [], 1);
sz = tmpEnd - tmpStart;
outBB = [tmpStart, sz];

end