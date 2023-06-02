function outPts = getPtsInsideView(inPts, view)
    idx = inPts(:,1) > view(2) & inPts(:,1) <= view(2) + view(4)...
        & inPts(:,2) > view(1) & inPts(:,2) <= view(1) + view(3);
    outPts = inPts(idx,:);
end