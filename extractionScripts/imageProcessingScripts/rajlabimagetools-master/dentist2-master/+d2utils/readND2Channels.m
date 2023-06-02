function channels = readND2Channels(fileND2)

    % Dictionary of commonly used channels and their names in Elements (in
    % Raj and Shaffer labs)
   channelMap = containers.Map({'Brightfield', 'DAPI', 'YFP', 'GFP', 'CY3', 'Cy3', 'A594', 'CY5', 'A647', '700', 'CY7','NIR'},...
                           {'trans'      , 'dapi', 'gfp', 'gfp', 'tmr', 'tmr', 'alexa', 'cy', 'cy'  , 'nir', 'nir','nir'});
%     channelMap = containers.Map({'DAPI','YFP','CY5','CY5_1','CY5_3','CY5_2',  'CY5_5','CY5_4','Brightfield',  'CY5a', 'CY5b', 'CY5c', 'CY5d','CY5e'},...
%                                 {'dapi','gfp','cy','cya',   'cyb',  'cyc',    'cyd',  'cye',  'trans',        'cya',  'cyb',  'cyc',  'cyd',  'cye'});
% warning('IanD modified readND2Channels.m')
    
    reader = bfGetReader(fileND2);
    omeMeta = reader.getMetadataStore();
    
    channels = cell(1, omeMeta.getPixelsSizeC(0).getValue());
    
    for i = 1:numel(channels)
        c = omeMeta.getChannelName(0, i-1);
        channels{i} = channelMap(c.toCharArray');
    end
    %channels = string(channels);
end