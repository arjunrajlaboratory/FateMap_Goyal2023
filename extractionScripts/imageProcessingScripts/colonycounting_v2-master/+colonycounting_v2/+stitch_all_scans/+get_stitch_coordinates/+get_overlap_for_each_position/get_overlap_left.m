function overlap_left = get_overlap_left(stitch_coords_center, stitch_coords_left, image_height)

    %%% NOTE: It is A LOT easier to understand this function if you refer
    %%% to the drawing included in the readme.

    %%% First, we want to get the coordinates of the corners of the
    %%% overlap. These coordinates will be in the reference frame of the
    %%% stitched image.

    % get the coords of the ul corner of the overlap:
    overlap_left.stitch.corner_ul_row = max(stitch_coords_center.corner_ul_row, stitch_coords_left.corner_ul_row);
    overlap_left.stitch.corner_ul_column = stitch_coords_center.corner_ul_column;

    % get the coords of the lr corner of the overlap:
    overlap_left.stitch.corner_lr_row = max(stitch_coords_center.corner_lr_row, stitch_coords_left.corner_lr_row);
    overlap_left.stitch.corner_lr_column = stitch_coords_left.corner_lr_column;

    % get the dimensions of the overlap:
    overlap_height = overlap_left.stitch.corner_lr_row - overlap_left.stitch.corner_ul_row;
    overlap_width = overlap_left.stitch.corner_lr_column - overlap_left.stitch.corner_ul_column;

    % get the coords of the ur corner of the overlap:
    [overlap_left.stitch.corner_ur_row, overlap_left.stitch.corner_ur_column] = get_ur_corner_coords(overlap_left.stitch.corner_ul_row, overlap_left.stitch.corner_ul_column, overlap_width);

    % get the coords of the ll corner of the overlap:
    [overlap_left.stitch.corner_ll_row, overlap_left.stitch.corner_ll_column] = get_ll_corner_coords(overlap_left.stitch.corner_ul_row, overlap_left.stitch.corner_ul_column, overlap_height);

    %%% Next, we want to convert the coordinates of the corners of the
    %%% overlap to the reference frame of an individual position.

    % if the center position is above the left position:
    if stitch_coords_center.corner_ul_row < stitch_coords_left.corner_ul_row

        % get coords of the ul corner of the overlap:
        overlap_left.position.corner_ul_row = image_height - overlap_height;

    % if the center position is below (or at the same height) as the left
    % position:
    elseif stitch_coords_center.corner_ul_row >= stitch_coords_left.corner_ul_row

        % get coords of the ul corner of the overlap:
        overlap_left.position.corner_ul_row = 1;

    end

    % get coords of the ul corner of the overlap:
    overlap_left.position.corner_ul_column = 1;

    % get coords of the lr corner of the overlap:
    [overlap_left.position.corner_lr_row, overlap_left.position.corner_lr_column] = get_lr_corner_coords(overlap_left.position.corner_ul_row, overlap_left.position.corner_ul_column, overlap_height, overlap_width);

    % get the coords of the ur corner of the overlap:
    [overlap_left.position.corner_ur_row, overlap_left.position.corner_ur_column] = get_ur_corner_coords(overlap_left.position.corner_ul_row, overlap_left.position.corner_ul_column, overlap_width);

    % get coords of the ll corner of the overlap:
    [overlap_left.position.corner_ll_row, overlap_left.position.corner_ll_column] = get_ll_corner_coords(overlap_left.position.corner_ul_row, overlap_left.position.corner_ul_column, overlap_height);

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