function centroids_keep = count_cells_in_scan(image, gaussian_sigma)

    %%% First, we want to smooth the image and find regional max. These
    %%% will be the centroids of the cells. 

    % set the gaussian filter to use:
    function_gaussian = @(block_struct) ...
        imgaussfilt(block_struct.data, gaussian_sigma);
    
    % smooth image with a gaussian (using block processing):
    image_blurred = blockproc(scale(image), [11000 11000], function_gaussian, 'BorderSize', [200 200]);

    % set the regional max to use:
    function_max = @(block_struct)...
        imregionalmax(block_struct.data);
    
    % sget regional max of the smoothed image (using block processing):
    image_max = blockproc(image_blurred, [11000 11000], function_max, 'BorderSize', [200 200]);

    % find connected components:
    CC = bwconncomp(image_max);
    
    % get centroids:
    centroids = regionprops(CC,'Centroid');
    centroids = [centroids.Centroid];
    centroids = reshape(centroids,2,[])'; 
    centroids = round(centroids); 
    
    %%% Next, we want to eliminate centroids where the image intensity is
    %%% below a threshold. This will get rid of any false positives. 

    % get the image intensity at every centroid:
    centroids_intensity = zeros(1, size(centroids, 1));
    for i = 1:size(centroids, 1)
        centroids_intensity(i) = image_blurred(centroids(i,2),centroids(i,1));
    end
    
    % create a histogram of the centroid intensities:
    figure('visible', 'off');
    histogram_handle = histogram(centroids_intensity);
    
    % get the x- and y-values of the histogram:
    histogram_y = histogram_handle.Values;
    histogram_x = histogram_handle.BinEdges;

    % format the x-values of the histogram to be from the center of each
    % bar:
    histogram_x = mean([histogram_x(1:end-1); histogram_x(2:end)]);
    
    % get indices where the slope is positive:
    histogram_y_derivative = diff(histogram_y);
    indices_slope_positive = find(histogram_y_derivative > 0) + 1; 
    
    % get indices past the max of the histogram:
    [~, index_max] = max(histogram_y);
    indices_beyond_max = find(histogram_x > histogram_x(index_max+1));
    
    % get the FIRST index where both conditions are true:
    index_threshold = indices_slope_positive(find(indices_slope_positive > indices_beyond_max(1), 1));
    
    % get the x-value at this index (this is the intensity threshold):
    threshold = histogram_x(index_threshold);

    % keep the centroids past the inflection point of the histogram:
    centroids_keep = centroids(centroids_intensity > threshold, :);

end