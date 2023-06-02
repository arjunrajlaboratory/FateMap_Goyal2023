function imOut = loadStitchD2(stitchFile, channel)
    tmpMat = matfile(stitchFile);
    if strcmp(channel, 'dapi')
        imOut = tmpMat.dapi;
    else
        fishLabels = tmpMat.fishLabels;
        idx = find(ismember(channel, fishLabels));
        tmpStitch = tmpMat.fishScans(1,idx);
        imOut = tmpStitch{:};
    end
end