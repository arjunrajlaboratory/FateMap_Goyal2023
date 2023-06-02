function boundaries = scale_coords_boundary_up(boundaries, scale_rows, scale_columns)
    
    % scale well boundaries:
    boundaries.well.stitch_original = scale_coords_boundary(boundaries.well.stitch_small, scale_rows, scale_columns);
    
    % scale colony boundaries:
    boundaries.colonies.stitch_original = scale_coords_boundary(boundaries.colonies.stitch_small, scale_rows, scale_columns);

end

function boundaries_original = scale_coords_boundary(boundaries, scale_rows, scale_columns)

    % get number of boundaries:
    num_boundaries = numel(boundaries);
    
    % add fields to store scaled coords:
    [boundaries_original(1:num_boundaries).coordinates_boundary] = deal([]);
    [boundaries_original(1:num_boundaries).coordinates_mask] = deal([]);

    % for each boundary:
    for i = 1:num_boundaries
       
        % get coords:
        coordinates_mask = boundaries(i).coordinates_mask;
        coordinates_boundary = boundaries(i).coordinates_boundary;
        
        % scale coords:
        coordinates_mask_original = colonycounting_v2.utilities.scale_coords(coordinates_mask, scale_rows, scale_columns);
        coordinates_boundary_original = colonycounting_v2.utilities.scale_coords(coordinates_boundary, scale_rows, scale_columns);
        
        % save:
        boundaries_original(i).coordinates_mask = coordinates_mask_original;
        boundaries_original(i).coordinates_boundary = coordinates_boundary_original;
        
    end
    
end