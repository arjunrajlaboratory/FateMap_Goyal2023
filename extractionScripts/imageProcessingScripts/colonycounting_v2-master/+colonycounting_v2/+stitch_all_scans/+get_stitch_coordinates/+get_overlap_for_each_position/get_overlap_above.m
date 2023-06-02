function overlap_above = get_overlap_above(stitch_coords_center, stitch_coords_above, image_width)

    %%% NOTE: It is A LOT easier to understand this function if you refer
    %%% to the drawing included in the readme.

    %%% First, we want to get the coordinates of the corners of the
    %%% overlap. These coordinates will be in the reference frame of the
    %%% stitched image.

    % get the coords of the ul corner of the overlap:
    overlap_above.stitch.corner_ul_row = stitch_coords_center.corner_ul_row;
    overlap_above.stitch.corner_ul_column = max(stitch_coords_center.corner_ul_column, stitch_coords_above.corner_ul_column);

    % get the coords of the lr corner of the overlap:
    overlap_above.stitch.corner_lr_row = stitch_coords_above.corner_lr_row;
    overlap_above.stitch.corner_lr_column = max(stitch_coords_center.corner_lr_column, stitch_coords_above.corner_lr_column);

    % get the dimensions of the overlap:
    overlap_height = overlap_above.stitch.corner_lr_row - overlap_above.stitch.corner_ul_row;
    overlap_width = overlap_above.stitch.corner_lr_column - overlap_above.stitch.corner_ul_column;

    % get the coords of the ur corner of the overlap:
    [overlap_above.stitch.corner_ur_row, overlap_above.stitch.corner_ur_column] = get_ur_corner_coords(overlap_above.stitch.corner_ul_row, overlap_above.stitch.corner_ul_column, overlap_width);

    % get the coords of the ll corner of the overlap:
    [overlap_above.stitch.corner_ll_row, overlap_above.stitch.corner_ll_column] = get_ll_corner_coords(overlap_above.stitch.corner_ul_row, overlap_above.stitch.corner_ul_column, overlap_height);

    %%% Next, we want to convert the coordinates of the corners of the
    %%% overlap to the reference frame of an individual position.

    % if the center position is to the left of the above position:
    if stitch_coords_center.corner_ul_column < stitch_coords_above.corner_ul_column

        % get coords of the ul corner of the overlap:
        overlap_above.position.corner_ul_column = image_width - overlap_width;

    % if the center position is to the right of (or at) the above position:
    elseif stitch_coords_center.corner_ul_column >= stitch_coords_above.corner_ul_column

        % get coords of the ul corner of the overlap:
        overlap_above.position.corner_ul_column = 1;

    end

    % get coords of the ul corner of the overlap:
    overlap_above.position.corner_ul_row = 1;

    % get coords of the lr corner of the overlap:
    [overlap_above.position.corner_lr_row, overlap_above.position.corner_lr_column] = get_lr_corner_coords(overlap_above.position.corner_ul_row, overlap_above.position.corner_ul_column, overlap_height, overlap_width);

    % get the coords of the ur corner of the overlap:
    [overlap_above.position.corner_ur_row, overlap_above.position.corner_ur_column] = get_ur_corner_coords(overlap_above.position.corner_ul_row, overlap_above.position.corner_ul_column, overlap_width);

    % get coords of the ll corner of the overlap:
    [overlap_above.position.corner_ll_row, overlap_above.position.corner_ll_column] = get_ll_corner_coords(overlap_above.position.corner_ul_row, overlap_above.position.corner_ul_column, overlap_height);

end

function [coords_ur_row, coords_ur_column] = get_ur_corner_coords(coords_ul_row, coords_ul_column, overlap_width)

    coords_ur_row = coords_ul_row;
    coords_ur_column = coords_ul_column + overlap_width;

end

function [coords_ll_row, coords_ll_column] = get_ll_corner_coords(coords_ul_row, coords_ul_column, overlap_height)

    coords_ll_row = coords_ul_row + overlap_height;
    coords_ll_column = coords_ul_column;

end

function [coords_lr_row, coords_lr_column] = get_lr_corner_coords(coords_ul_row, coords_ul_column, overlap_height, overlap_width)

    coords_lr_row = coords_ul_row + overlap_height;
    coords_lr_column = coords_ul_column + overlap_width;

end