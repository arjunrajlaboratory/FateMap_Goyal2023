function [] = findBarcodeProbesHCR(path)
%%Assign fasta file names to array
    %path = '/Volumes/B_MR/timeMachineNGS/20190506_TM22_run2';
    fastaFiles = dir(fullfile(path, 'fastaFiles/barcode*.fa'));
    barcodeFileNames = {fastaFiles.name};
    barcodeFilePath = {fastaFiles.folder};
    barcodeFastaFiles = join([barcodeFilePath; barcodeFileNames], '/', 1);
    [folder, baseFileName, extension] = cellfun(@fileparts, barcodeFastaFiles, 'UniformOutput', false);
    outFiles = strcat(fullfile(path, 'probeDesignHCR', baseFileName), {'_HCR_42mers'});
    pseudogeneMaskOption = true;
    %% Make probeDesign subdirectory if it doesn't already exist
    if ~exist(fullfile(path, 'probeDesignHCR'), 'dir')
        mkdir(fullfile(path, 'probeDesignHCR'))
    end
    %%
    for i = 1:length(barcodeFastaFiles)
        findprobesHD(barcodeFastaFiles{i}, 2, 'outfilename', outFiles{i}, 'oligolength', 42, 'spacerlength', 0, 'allowableGibbsFE', [-65,-45], 'targetGibbsFE', -55, 'pseudogenemask', pseudogeneMaskOption, 'masksequences', '/Volumes/B_MR2/timeMachineNGS/vectorMaps/lentiefs_gfp_m20_0mismatch_a2.fasta', 'nummatchesformask', 25)
    end
    %% old script below
%     for i = 1:length(barcodeFastaFiles)
%         err_count = 0
%         while err_count < 3
%             try
%                 findprobesHD(barcodeFastaFiles{i}, 2, 'outfilename', outFiles{i}, 'oligolength', 42, 'spacerlength', 0, 'allowableGibbsFE', [-65,-45], 'targetGibbsFE', -55, 'pseudogenemask', pseudogeneMaskOption, 'masksequences', '/Volumes/B_MR2/timeMachineNGS/vectorMaps/lentiefs_gfp_m20_0mismatch_a2.fasta', 'nummatchesformask', 25)
%                 err_count = 3
%             catch
%                 err_count = err_count + 1
%             end
%         end
%     end
end