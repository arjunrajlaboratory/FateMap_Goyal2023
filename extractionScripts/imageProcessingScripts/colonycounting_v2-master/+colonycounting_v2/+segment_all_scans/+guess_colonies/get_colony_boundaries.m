function boundaries = get_colony_boundaries(boundaries, cells_x, cells_y, colony_assignment, stitch_small)

    % get number of colonies:
    num_colonies = numel(unique(colony_assignment));
    
    % for each colony:
    for i = 1:num_colonies
        
        % get indices of cells in the colony:
        temp_indices = find(colony_assignment == i);

        % get cells in the colony:
        temp_cells_colony_x = cells_x(temp_indices);
        temp_cells_colony_y = cells_y(temp_indices);

        % if there are more than 10 cells in the colony:
        if numel(temp_cells_colony_x) > 10

            % get the convex hull of the colony:
            convex_hull = convhull(temp_cells_colony_x, temp_cells_colony_y);

            % get the boundaries of the convex hull:
            temp_colony_boundary_x = temp_cells_colony_x(convex_hull);
            temp_colony_boundary_y = temp_cells_colony_y(convex_hull);
            
            % get a mask of the colony:
            temp_mask = poly2mask(temp_colony_boundary_x, temp_colony_boundary_y, size(stitch_small, 1), size(stitch_small, 2));
            
            % get coordinates of the mask:
            [temp_colony_mask_x, temp_colony_mask_y] = find(temp_mask == 1);

            % save:
            temp.status = 'keep';
            temp.coordinates_boundary = [temp_colony_boundary_x, temp_colony_boundary_y];
            temp.coordinates_mask = [temp_colony_mask_x, temp_colony_mask_y];
            boundaries = colonycounting_v2.utilities.add_entry_to_structure(temp, boundaries);

        end
        
    end

end