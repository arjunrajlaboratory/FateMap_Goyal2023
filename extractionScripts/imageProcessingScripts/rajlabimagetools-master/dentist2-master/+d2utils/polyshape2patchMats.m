function [xmat, ymat] = polyshape2patchMats(polyIn)
    tmp = {polyIn.Vertices};
    xCell = cellfun(@(x) x(:,1), tmp, 'UniformOutput', false);
    yCell = cellfun(@(x) x(:,2), tmp, 'UniformOutput', false);
    maxPoints = max(cellfun(@height, xCell));
    xCell2 = cellfun(@(x) [x; repmat(x(1), maxPoints-numel(x), 1)], xCell, 'UniformOutput', false);
    yCell2 = cellfun(@(x) [x; repmat(x(1), maxPoints-numel(x), 1)], yCell, 'UniformOutput', false);
    xmat = [xCell2{:}];
    ymat = [yCell2{:}];
end