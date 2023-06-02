function scan = get_shift_for_each_position(scan)

    % for each position:
    for i = 1:numel(scan.stitch_coords)

        % determine the shift (in x and y) due to the row:
        shift_in_x_row = scan.stitch_coords(i).position_num_column * scan.shift_row.column;
        shift_in_y_row = scan.stitch_coords(i).position_num_column * scan.shift_row.row;

        % determine the shift (in x and y) due to the column:
        shift_in_x_column = scan.stitch_coords(i).position_num_row * scan.shift_column.column;
        shift_in_y_column = scan.stitch_coords(i).position_num_row * scan.shift_column.row;

        % determine the total shift:
        scan.stitch_coords(i).shift_in_x = shift_in_x_row + shift_in_x_column;
        scan.stitch_coords(i).shift_in_y = shift_in_y_row + shift_in_y_column;

    end

end