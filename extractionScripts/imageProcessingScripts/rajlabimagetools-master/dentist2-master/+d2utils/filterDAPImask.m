function [outMask, rp]= filterDAPImask(inMask, minObjectSize)
    CC = bwconncomp(inMask);
    rp = regionprops(CC);
    area = [rp.Area];

    idx = area < minObjectSize;
    rp = rp(idx);

    filterList = CC.PixelIdxList(idx);
    outMask = inMask;
    for i = 1:numel(filterList)
        outMask(filterList{i}) = 0;
    end
end