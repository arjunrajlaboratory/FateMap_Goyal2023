function tmpImage = ND2tileReader(inFile, varargin)
    p = inputParser;

    p.addRequired('inFile', @ischar);
    p.addParameter('tile', 1, @isnumeric);
    
    p.addParameter('channel', 1, @isnumeric);
    p.addParameter('Z', 1, @isnumeric);
    
    p.parse(inFile, varargin{:});
    
    reader = bfGetReader(p.Results.inFile);
    
    reader.setSeries(p.Results.tile-1);
    iPlane = reader.getIndex(p.Results.Z - 1, p.Results.channel - 1, 0) + 1;
    
    tmpImage  = bfGetPlane(reader, iPlane);
end
