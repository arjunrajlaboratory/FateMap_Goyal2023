function cells_and_boundaries = get_cells_within_boundaries(cells, boundaries)

    % get number of boundaries:
    num_boundaries = numel(boundaries.stitch_original);
    
    % create structure to store cells and boundaries:
    [cells_and_boundaries(1:num_boundaries).boundary_stitch] = deal([]);
    [cells_and_boundaries(1:num_boundaries).boundary_stitch_small] = deal([]);
    [cells_and_boundaries(1:num_boundaries).cells_stitch] = deal([]);
    [cells_and_boundaries(1:num_boundaries).cells_stitch_small] = deal([]);
    
    % for each boundary:
    for i = 1:num_boundaries

        % get boundary coordinates:
        boundary_row = boundaries.stitch_original(i).coordinates_boundary(:,1);
        boundary_col = boundaries.stitch_original(i).coordinates_boundary(:,2);
        
        % get cell coordinates:
        cells_row = cells.stitch(:,1);
        cells_col = cells.stitch(:,2);
        
        % get row indices of cells that fall within the boundary:
        rows = inpolygon(cells_col, cells_row, boundary_col, boundary_row);

        % save centroids that overlap with object:
        cells_and_boundaries(i).boundary_stitch = boundaries.stitch_original(i).coordinates_boundary;
        cells_and_boundaries(i).boundary_stitch_small = boundaries.stitch_small(i).coordinates_boundary;
        cells_and_boundaries(i).cells_stitch = cells.stitch(rows, :);
        cells_and_boundaries(i).cells_stitch_small = cells.stitch_small(rows, :);
        
    end

end