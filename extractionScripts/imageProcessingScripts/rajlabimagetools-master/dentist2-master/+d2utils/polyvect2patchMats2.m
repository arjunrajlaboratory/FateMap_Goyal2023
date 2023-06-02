function [xmat, ymat, colorMat] = polyvect2patchMats2(polyIn, colors)
    singleRegionIdx = [polyIn.NumRegions] == 1;
    tmp = {polyIn.Vertices};
    maxPoints = max(cellfun(@height, tmp));
    tmp = tmp(singleRegionIdx);
    tmp = cellfun(@(x) [x; repmat(x(1,:), maxPoints-height(x), 1)], tmp, 'UniformOutput', false);
    xmat = cell2mat(cellfun(@(x) x(:,1), tmp, 'UniformOutput', false));
    ymat = cell2mat(cellfun(@(x) x(:,2), tmp, 'UniformOutput', false));
    colorMat = colors(singleRegionIdx, :);
    if any(~singleRegionIdx)
        tmp = polyIn(~singleRegionIdx);
        tmpColors = colors(~singleRegionIdx, :);
        for i = 1:numel(tmp)
            tmpPoly = regions(tmp(i));
            tmpVertices = {tmpPoly.Vertices};
            tmpVertices = cellfun(@(x) [x; repmat(x(1,:), maxPoints-height(x), 1)], tmpVertices, 'UniformOutput', false);
            tmpXmat = cell2mat(cellfun(@(x) x(:,1), tmpVertices, 'UniformOutput', false));
            tmpYmat = cell2mat(cellfun(@(x) x(:,2), tmpVertices, 'UniformOutput', false));
            xmat = [xmat, tmpXmat];
            ymat = [ymat, tmpYmat];
            colorMat = [colorMat; repmat(tmpColors(i, :), width(tmpXmat),1)];
        end
    end
end