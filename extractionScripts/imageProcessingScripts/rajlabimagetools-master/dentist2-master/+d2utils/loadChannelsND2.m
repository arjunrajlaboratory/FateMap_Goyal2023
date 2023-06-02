function channelsOut = loadChannelsND2(infile, channelsIn)
    if isempty(channelsIn)
        channelsOut = d2utils.readND2Channels(infile);
    elseif numel(channelsIn) ~= numel(unique(channelsIn)) %Consider replacing with function to auto rename duplicates
        error('Please specify only unique channel names.')
    else
        channelsOut = channelsIn;
    end
end