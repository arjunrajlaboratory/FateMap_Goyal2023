function centroids_keep = count_cells_in_scan_adaptThresh(image, sensitivity, minNucleusSize, scaleOption)

    %%% Need to update script to incorporate scaling option 
 
    % binarize image to find centroids:
%     if scaleOption == true
%         function_binarize = @(block_struct)...
%             imbinarize(scale(block_struct.data), adaptthresh(scale(block_struct.data), sensitivity, 'ForegroundPolarity','bright')); 
%     else
%         function_binarize = @(block_struct)...
%             imbinarize(block_struct.data, adaptthresh(block_struct.data, sensitivity, 'ForegroundPolarity','bright'));
%     end
    
    % binarize image to find centroids:
    function_binarize = @(block_struct)...
        imbinarize(scale(block_struct.data), adaptthresh(scale(block_struct.data), sensitivity, 'ForegroundPolarity','bright')); 
    
    % find centroids (using block processing):
    image_binarize = blockproc(image, [11000 11000], function_binarize, 'BorderSize', [200 200]);

    % find connected components and remove centroid smaller than threshold:
    CC = bwconncomp(image_binarize, 4);
        
    rp = regionprops(CC);
    area = [rp.Area];
    centroids = [rp.Centroid];
    centroids = reshape(centroids,2,[])'; 
    centroids = round(centroids); 

    idx = area > minNucleusSize; % Get rid of small stuff

    centroids_keep = centroids(idx, 1:end);
 
end