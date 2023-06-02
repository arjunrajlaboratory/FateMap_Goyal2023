function [xmat, ymat] = masktable2patchMats(maskTableIn)
    xCell = splitapply(@(x){x}, maskTableIn.x, findgroups(maskTableIn.maskID));
    yCell = splitapply(@(x){x}, maskTableIn.y, findgroups(maskTableIn.maskID));
    maxPoints = max(cellfun(@height, xCell));
    xCell2 = cellfun(@(x) [x; repmat(x(1), maxPoints-numel(x), 1)], xCell', 'UniformOutput', false);
    yCell2 = cellfun(@(x) [x; repmat(x(1), maxPoints-numel(x), 1)], yCell', 'UniformOutput', false);
    xmat = [xCell2{:}];
    ymat = [yCell2{:}];
end