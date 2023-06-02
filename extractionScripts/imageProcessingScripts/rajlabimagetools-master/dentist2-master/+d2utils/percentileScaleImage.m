function outIm = percentileScaleImage(img, percentiles, scaleFactors)
    tmpIm = im2single(img);
    plower = mean(prctile(tmpIm, percentiles(1)));
    pupper = mean(prctile(tmpIm, percentiles(2)));
    maxIntensity = min(pupper*scaleFactors(2), 1);
    outIm = tmpIm - plower*scaleFactors(1);
    outIm = outIm ./ maxIntensity;
    
    outIm(outIm<0) = 0;
    outIm(outIm>1) = 1;
end