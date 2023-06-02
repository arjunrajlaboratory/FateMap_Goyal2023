function [stitch_cropped, cells_algorithm1, cells_algorithm2] = get_example(stitch)
    
    % set the default parameter for orgiginal algorithm:
    gaussian_sigma = 7;
    
    %Set default parameters for adaptiveThreshold algorithm
    sensitivity = 0.3;
    minNucleusSize = 100;
    scaleOption = true;
    
    % set the default position for cropping image:
    position_original.x = round(size(stitch, 2) / 2);
    position_original.y = round(size(stitch, 1) / 2);
    position_original.width = 2000;
    position_original.height = 2000;
    
    % crop the original stitch:
    stitch_cropped = imcrop(stitch, [position_original.x, position_original.y, position_original.width, position_original.height]);
    
    %Run the original algorithm for identifying cells:
    cells_algorithm1 = colonycounting_v2.count_cells_all_scans.count_cells_in_scan(stitch_cropped, gaussian_sigma);
    
    %Run the adaptiveThreshold algorithm for identifying cells:
    cells_algorithm2 = colonycounting_v2.count_cells_all_scans.count_cells_in_scan_adaptThresh(stitch_cropped, sensitivity, minNucleusSize, scaleOption);
    
end