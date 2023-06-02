function [image_main, image_below, image_right] = get_images_from_tile_numbers(tile_row, tile_column, matrix_of_positions, list_images, path_images)

    % get the position of the image:
    position_num_main = matrix_of_positions(tile_row, tile_column);
    
    % get the position of the image below:
    position_num_below = matrix_of_positions(tile_row + 1, tile_column);
    
    % get the position of the image right:
    position_num_right = matrix_of_positions(tile_row, tile_column + 1);
    
    % convert position numbers to strings:
    position_num_main_string = ['s' num2str(position_num_main, '%.0f')];
    position_num_below_string = ['s' num2str(position_num_below, '%.0f')];
    position_num_right_string = ['s' num2str(position_num_right, '%.0f')];
    
    % get name and path of images:
    image_info_main = colonycounting_v2.utilities.get_structure_results_containing_string(list_images, 'name', position_num_main_string);
    image_info_below = colonycounting_v2.utilities.get_structure_results_containing_string(list_images, 'name', position_num_below_string);
    image_info_right = colonycounting_v2.utilities.get_structure_results_containing_string(list_images, 'name', position_num_right_string);
    
    % load the images:
    image_main = readmm(fullfile(path_images, image_info_main(1).name));
    image_main = image_main.imagedata;
    image_below = readmm(fullfile(path_images, image_info_below(1).name));
    image_below = image_below.imagedata;
    image_right = readmm(fullfile(path_images, image_info_right(1).name));
    image_right = image_right.imagedata;

end