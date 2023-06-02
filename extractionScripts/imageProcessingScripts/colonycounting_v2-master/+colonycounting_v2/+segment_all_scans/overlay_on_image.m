function image_annotated = overlay_on_image(image, cells)

    % create rgb-version of image:
    image_annotated = repmat(image, 1, 1, 3);
    
    % set marker size:
    marker_size = 1;
    
    %%% NOTE: I add the various cells and boundaries in a particular order
    %%% so that things are most visible. 
    
    %%% Next, we want to plot ALL cells.
    
    % set color:
    color.cells.all = [1 1 1];
    
    % add cells to image:
    image_annotated = colonycounting_v2.utilities.add_cells_to_image(image_annotated, cells.all.stitch_small, color.cells.all, marker_size);
    
    %%% Next, we want to plot the well boundary.
    
    % set color for the well boundary:
    color.well = [0 0 1];
    
    % for each well boundary:
    for i = 1:numel(cells.well)
        
        % add boundary to image:
        image_annotated = colonycounting_v2.utilities.add_boundary_to_image(image_annotated, cells.well(i).boundary_stitch_small, color.well, marker_size);
        
    end
    
    %%% Next, we want to plot all cells within the well.
    
    % for each well:
    for i = 1:numel(cells.well)
        
        % add cells to image:
        image_annotated = colonycounting_v2.utilities.add_cells_to_image(image_annotated, cells.well(i).cells_stitch_small, color.well, marker_size);
    
    end
    
    %%% Next, we want to plot all colony boundaries:
    
    % get number of colonies:
    num_colonies = numel(cells.colonies);
    
    % set colors for each colony:
    color.colonies = distinguishable_colors(num_colonies, {'w', 'k', 'b'});
    
    % for each colony:
    for i = 1:num_colonies
        % Skip colonies with fewer than 6 points on its perimeter
        % (likely bad segments)
        
        if size(cells.colonies(i).boundary_stitch_small, 1) > 5
            % add boundary to image:
            image_annotated = colonycounting_v2.utilities.add_boundary_to_image(image_annotated, cells.colonies(i).boundary_stitch_small, color.colonies(i, :), marker_size);
        end   
    end
    
    %%% Next, we want to plot all cells assigned to colonies:
    
    % for each colony:
    for i = 1:num_colonies
        % Skip colonies with fewer than 6 points on its perimeter
        % (likely bad segments)
        
        if size(cells.colonies(i).boundary_stitch_small, 1) > 5
            % add cells to image:
            image_annotated = colonycounting_v2.utilities.add_cells_to_image(image_annotated, cells.colonies(i).cells_stitch_small, color.colonies(i, :), marker_size);
        end
    end

end