function [num_tiles_row, num_tiles_column] = enter_scan_dimensions(path_folder, name_scan)

    % ask user for the scan dimensions:
    prompts = {path_folder, 'Number of Rows', 'Number of Columns'};
    defaults = {name_scan, '', ''};
    answer = inputdlg(prompts, 'Enter Tile Dimensions', [1 200], defaults);
    num_tiles_row = str2double(answer{2});
    num_tiles_column = str2double(answer{3});

end