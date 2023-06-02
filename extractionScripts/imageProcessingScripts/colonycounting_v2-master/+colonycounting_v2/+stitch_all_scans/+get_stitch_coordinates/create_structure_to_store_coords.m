function stitch_coords = create_structure_to_store_coords(scan)

    % get number of positions:
    num_positions = scan.num_tiles_row * scan.num_tiles_column;

    % create array with list of position numbers:
    array_of_positions = 1:num_positions;

    % arrange position numbers into matrix (same order that scope acquires images):
    matrix_of_positions = reshape(array_of_positions, scan.num_tiles_row, scan.num_tiles_column);

    % create structure:
    [stitch_coords(1:num_positions).position_num_linear] = deal(0);
    [stitch_coords(1:num_positions).position_num_row] = deal(0);
    [stitch_coords(1:num_positions).position_num_column] = deal(0);
    [stitch_coords(1:num_positions).shift_in_x] = deal(0);
    [stitch_coords(1:num_positions).shift_in_y] = deal(0);
    [stitch_coords(1:num_positions).corner_ul_row] = deal(0);
    [stitch_coords(1:num_positions).corner_ul_column] = deal(0);
    [stitch_coords(1:num_positions).corner_lr_row] = deal(0);
    [stitch_coords(1:num_positions).corner_lr_column] = deal(0);
    [stitch_coords(1:num_positions).overlap_left] = deal(0);
    [stitch_coords(1:num_positions).overlap_above] = deal(0);

    % fill in structure:
    for i = 1:num_positions

        % save the linear position number (same number scope uses to label
        % positions):
        stitch_coords(i).position_num_linear = i;

        % save the row and columb numbers:
        [stitch_coords(i).position_num_row, stitch_coords(i).position_num_column] = find(matrix_of_positions == stitch_coords(i).position_num_linear);

    end

end