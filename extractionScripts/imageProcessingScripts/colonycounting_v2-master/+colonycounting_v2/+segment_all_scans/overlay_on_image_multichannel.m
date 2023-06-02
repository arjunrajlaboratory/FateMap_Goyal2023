function image_annotated = overlay_on_image_multichannel(image, cells, groups)

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
    
    %%% Next, we want to plot all colony boundaries colored by group as well as the cells within boundaries:
    
    % set colors for each group:
    color.colonies = distinguishable_colors(numel(groups), {'w', 'k', 'b'});
   
 
    for i = 1:numel(groups)
        % get number of colonies:
        num_colonies = numel(cells.(groups{i}).colonies);
            % for each colony:
        for k = 1:num_colonies
        
            % add boundary to image:
            image_annotated = colonycounting_v2.utilities.add_boundary_to_image(image_annotated, cells.(groups{i}).colonies(k).boundary_stitch_small, color.colonies(i, :), marker_size);
                  
            % add cells to image:
            image_annotated = colonycounting_v2.utilities.add_cells_to_image(image_annotated, cells.(groups{i}).colonies(k).cells_stitch_small, color.colonies(i, :), marker_size);
   
        end
       
    end

end