function scan = get_overlap_for_each_position(scan)

    % for each position:
    for i = 1:numel(scan.stitch_coords)
        
        % get the position numbers:
        position_num_column = scan.stitch_coords(i).position_num_column;
        position_num_row = scan.stitch_coords(i).position_num_row;

        % get the stitch coords for the position:
        stitch_coords_center = colonycounting_v2.utilities.get_structure_results_matching_number(scan.stitch_coords, 'position_num_column', position_num_column);
        stitch_coords_center = colonycounting_v2.utilities.get_structure_results_matching_number(stitch_coords_center, 'position_num_row', position_num_row);

        %%% Note that we look for neighbors to the left and top because of
        %%% the order that images are added to the scan. Images are added
        %%% starting with the position in the lower right corner of the
        %%% stitch. Subseqeunt images are placed on top. 
        
        % if there is a position to the left:
        if position_num_column > 1
           
            % get the stitch coords for the position to the left:
            stitch_coords_left = colonycounting_v2.utilities.get_structure_results_matching_number(scan.stitch_coords, 'position_num_column', position_num_column - 1);
            stitch_coords_left = colonycounting_v2.utilities.get_structure_results_matching_number(stitch_coords_left, 'position_num_row', position_num_row);
            
            % get overlap to the left:
            overlap_left = colonycounting_v2.stitch_all_scans.get_stitch_coordinates.get_overlap_for_each_position.get_overlap_left(stitch_coords_center, stitch_coords_left, scan.image_height);
          
        % otherwise:
        else
            
            % set left overlap to none:
            overlap_left.stitch = 'none';
            overlap_left.position = 'none';
            
        end
        
        % if there is a position above:
        if position_num_row > 1
            
            % get the stitch coords for the position above:
            stitch_coords_above = colonycounting_v2.utilities.get_structure_results_matching_number(scan.stitch_coords, 'position_num_column', position_num_column);
            stitch_coords_above = colonycounting_v2.utilities.get_structure_results_matching_number(stitch_coords_above, 'position_num_row', position_num_row - 1);
            
            % create structure to store overlap:
            overlap_above = colonycounting_v2.stitch_all_scans.get_stitch_coordinates.get_overlap_for_each_position.get_overlap_above(stitch_coords_center, stitch_coords_above, scan.image_width);
            
        % otherwise:
        else
            
            % set the above overlap to none:
            overlap_above.stitch = 'none';
            overlap_above.position = 'none';
            
        end
        
        % save overlap:
        scan.stitch_coords(i).overlap_left = overlap_left;
        scan.stitch_coords(i).overlap_above = overlap_above;
        
    end

end