function channelsOut = loadChannelsTiff(infile, channelsIn)
    stackInfo = imfinfo(infile);
    if isempty(channelsIn)
        error('To load multichannel tiff file, please specify channels with the "channels" parameter.')
    elseif numel(channelsIn) ~= numel(unique(channelsIn)) %Consider replacing with function to auto rename duplicates
        error('Input channel names contain duplicates. Please specify only unique channel names.')
    elseif numel(stackInfo) ~= numel(channelsIn)
        error('Pleaes input 1 channel name for each channel in your scan file.\nThe input scan file: %s contains %d channels.\n%d channel names were specified\n', infile, numel(stackInfo), numel(channelsIn))
    else
        channelsOut = channelsIn;
    end
end