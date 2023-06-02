function outIm = percentileAdjustImage(img, percentiles, scaleFactors)
    tmpIm = im2single(img);
    plower = mean(prctile(tmpIm, percentiles(1)));
    pupper = mean(prctile(tmpIm, percentiles(2)));
    outIm = imadjust(tmpIm, [plower*scaleFactors(1) min(pupper*scaleFactors(2), 1)]);
    outIm = im2uint16(outIm);
end