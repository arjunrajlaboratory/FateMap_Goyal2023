function [] = nd2ScanToTiff(varargin)
% Script to take scan acquired using Nikon elements and convert to
% tiff format compatible with Rajlab image tools. Require bfmatlab code is in
% your matlab path for reading nd2 files. The bfmatlab code can be
% downloaded here https://docs.openmicroscopy.org/bio-formats/5.7.3/developers/matlab-dev.html

% Example of usages:
% nd2ScanToTiff('inDir', 'eduardoTestScan/plateA')
% nd2ScanToTiff('inDir', 'eduardoTestScan/plateA', 'outdir', 'eduardoTestScan/plateA/splitPoints')
p = inputParser;

p.addParameter('inDir', '', @ischar);
p.addParameter('outDir', '', @ischar);

p.parse(varargin{:});

if ~isempty(p.Results.inDir)
    inDir = p.Results.inDir;
else
    inDir = pwd;
end

if ~isempty(p.Results.outDir)
    outDir = p.Results.outDir;
else
    outDir = fullfile(inDir, 'renamed');
end

% Make directory for tiff images
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

inFiles = dir(fullfile(inDir, '*.nd2'));
inFiles = {inFiles.name};
nScans = length(inFiles);

%Check if any scan files already exit
if isempty(dir(fullfile(outDir, 'Scan*.TIF')))
    scanCount = 0;
else
    currentScans = dir(fullfile(outDir, 'Scan*.TIF'));
    currentScans = {currentScans.name};
    currentScanCount = regexp(currentScans, '(?<=Scan)\d{3}?', 'match');
    currentScanCount = vertcat(currentScanCount{:});
    currentScanCount = str2num(cell2mat(currentScanCount));
    scanCount = 0 + max(currentScanCount); 
end

%Loop through ndFiles, saving each as a seperate Scan00#
for i = 1:nScans
    reader = bfGetReader(fullfile(inDir, inFiles{i}));
    omeMeta = reader.getMetadataStore();
    
    % Get metadata from .nd2 file
    stackSizeC = omeMeta.getPixelsSizeC(0).getValue(); % number of wavelength channels
    imageN = omeMeta.getImageCount();              % number of positions
%   stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices

    % Arrays if stacks are being acquired at each position
%   wave_numb = repmat(1:stackSizeC,1,stackSizeZ);     % repeating vector of channel #
%   Z_numb    = repmat(1:stackSizeZ,1,stackSizeC);     % repeating vector of Z #
    
    % For this scan, see if any TIF files exist in output directory
    if isempty(dir(sprintf('Scan%03d_*.TIF', i)))
         imageCount = 1;
    else
        currentFiles = dir(sprintf('Scan%03d_*.TIF', scanCount));
        currentFiles = {currentFiles.name};
        currentFileCount = regexp(currentFiles, '(?<=s)\d+', 'match');
        currentFileCount = vertcat(currentFileCount{:});
        currentFileCount = str2num(cell2mat(currentFileCount));
        imageCount = 1 + max(currentFileCount); 
    end
    
    for ii = imageCount:imageN % Loop through each image
        reader.setSeries(ii-1)

        for iii = 1:stackSizeC  % (# channels)
            % Read plane from series iSeries at Z, C, T coordinates (iZ, iC, iT)
  %         iPlane = reader.getIndex(Z_numb(iv)-1, wave_numb(iii) - 1, 0) + 1; 
            
            iPlane = reader.getIndex(0, iii - 1, 0) + 1; 
            stack_fig  = bfGetPlane(reader, iPlane);

            outputFileName = fullfile(outDir, sprintf('Scan%03d_w%d_s%d_t1.TIF', i+scanCount, iii, ii));
            imwrite(stack_fig, outputFileName)

        end
    end
    sprintf('Done with file %s', fullfile(inDir, inFiles{i}))
end
sprintf('Done with directory %s', inDir)
end

