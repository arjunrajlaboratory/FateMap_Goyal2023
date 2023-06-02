function stitch_info = get_stitch_coordinates(stitch_info)

    % for each folder:
    for i = 1:numel(stitch_info)
       
        % for each scan:
        for j = 1:numel(stitch_info(i).scan_info)
            
            %%% First, we want to set up a structure to store all the coordinates
            %%% and other information necessary to stitching the scan.
            
            stitch_info(i).scan_info(j).stitch_coords = colonycounting_v2.stitch_all_scans.get_stitch_coordinates.create_structure_to_store_coords(stitch_info(i).scan_info(j));

            %%% Next, we want to determine by how much to shift each pixel in x
            %%% and y (accounting for both the shift due to the row and the shift
            %%% due to the column). 
            
            stitch_info(i).scan_info(j) = colonycounting_v2.stitch_all_scans.get_stitch_coordinates.get_shift_for_each_position(stitch_info(i).scan_info(j));
            
            %%% Next, we want to convert the shift amounts into coordinates in the
            %%% reference frame of the stitched image. The easiest way to do this
            %%% is to offset the shifts so that they start at 1. 

            stitch_info(i).scan_info(j) = colonycounting_v2.stitch_all_scans.get_stitch_coordinates.get_coords_for_each_position(stitch_info(i).scan_info(j));
        
            %%% Next, we want to determine, for each position, the parts of the
            %%% image that overlap with another position. This will be used later
            %%% on to make sure cells are not double-counted. 
            
            stitch_info(i).scan_info(j) = colonycounting_v2.stitch_all_scans.get_stitch_coordinates.get_overlap_for_each_position(stitch_info(i).scan_info(j));
            
        end
        
    end

end