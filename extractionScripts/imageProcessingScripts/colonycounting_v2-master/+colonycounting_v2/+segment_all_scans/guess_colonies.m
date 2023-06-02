function boundaries = guess_colonies(cells, boundaries, stitch_small, boundary_well)

    % exclude cells outside the well:
    cells_inside = inpolygon(cells(:,1), cells(:,2), boundary_well.coordinates_boundary(:,1), boundary_well.coordinates_boundary(:,2));

    % get cell coordinates:
    cells_x = cells(cells_inside,1);
    cells_y = cells(cells_inside,2);
    
    % get cell coordinates rounded:
    cells_x_rounded = max(min(round(cells_x), size(stitch_small, 2)), 0);
    cells_y_rounded = max(min(round(cells_y), size(stitch_small, 1)), 0);

    % get number of cells:
    num_cells = numel(cells_x);
    
    % create binary image with cell positions:
    image = zeros(size(stitch_small));
    for i = 1:num_cells
        y = cells_y_rounded(i);
        x = cells_x_rounded(i);
        image(y, x) = image(y, x) + 1; 
    end

    % get density image:
    image_density = colonycounting_v2.segment_all_scans.guess_colonies.create_density_image(image, 5);

    % thresold density image:
    image_density_centroids = image_density;
    image_density_centroids(image_density_centroids < 20) = 0;
    
    % threshold the image using the density image (ignore cells in regions
    % of low density):
    image_threshold = image;
    image_threshold(~(image_density_centroids & image_threshold)) = 0;
    
    % get the centroids of the colonies:
    [colony_centroids_x, colony_centroids_y] = colonycounting_v2.segment_all_scans.guess_colonies.get_colony_centroids(image_threshold);
    
    % get number of colonies:
    num_colonies = numel(colony_centroids_x);
    
    % assign cells to colonies:
    colony_assignment = colonycounting_v2.segment_all_scans.guess_colonies.assign_cells_to_centroids(cells_x_rounded, cells_y_rounded, image_threshold, colony_centroids_x, colony_centroids_y);
    
    % use the cell assignments to get colony boundaries:
    boundaries = colonycounting_v2.segment_all_scans.guess_colonies.get_colony_boundaries(boundaries, cells_x_rounded, cells_y_rounded, colony_assignment, stitch_small);
    
end