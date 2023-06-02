function [outputOligos] = findBarcodeProbesClamp(path)
%%Assign fasta file names to array
fastaFiles = dir(fullfile(path, 'fastaFiles/barcode*.fa'));
barcodeFileNames = {fastaFiles.name};
barcodeFilePath = {fastaFiles.folder};
barcodeFastaFiles = join([barcodeFilePath; barcodeFileNames], '/', 1);
[folder, baseFileName, extension] = cellfun(@fileparts, barcodeFastaFiles, 'UniformOutput', false);
outFiles = strcat(fullfile(path, 'probeDesign', baseFileName), {'_clampFISH_backbones'});
%% Make probeDesign subdirectory if it doesn't already exist
if ~exist(fullfile(path, 'probeDesign'))
    mkdir(fullfile(path, 'probeDesign'))
end
%% Call findprobesHD on each file in fastaFiles array.
for i = 0:length(barcodeFastaFiles)
    err_count = 0
    while err_count < 3
        try
            findprobesHD(barcodeFastaFiles{i}, 3, 'outfilename', outFiles{i}, 'oligolength', 30, 'spacerlength', 0, 'allowableGibbsFE', [-50,-30], 'targetGibbsFE', -40, 'masksequences', '/Volumes/B_MR/timeMachineNGS/vectorMaps/lentiefs_gfp_m20_0mismatch_a2.fasta', 'nummatchesformask', 20)
            err_count = 3
        catch
            err_count = err_count + 1
        end
    end
end
