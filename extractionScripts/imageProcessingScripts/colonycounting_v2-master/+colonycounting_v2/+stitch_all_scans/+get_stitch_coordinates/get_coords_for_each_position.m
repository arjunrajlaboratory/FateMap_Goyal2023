function scan = get_coords_for_each_position(scan)

    % load an image to get the image size:
    temp_image = imread(fullfile(scan.path_folder, scan.images(1).list_images(1).name));
    image_height = size(temp_image, 1);
    image_width = size(temp_image, 2);

    % get minimum shift amounts:
    min_shift_x = extractfield(scan, 'stitch_coords');
    min_shift_x = min(extractfield(min_shift_x{1}, 'shift_in_x'));
    min_shift_y = extractfield(scan, 'stitch_coords');
    min_shift_y = min(extractfield(min_shift_y{1}, 'shift_in_y'));

    % for each position:
    for i = 1:numel(scan.stitch_coords)

        % offset the shift amounts to create the coordinates of the upperleft
        % corner of each image:
        scan.stitch_coords(i).corner_ul_row = scan.stitch_coords(i).shift_in_y - min_shift_y + 1;
        scan.stitch_coords(i).corner_ul_column = scan.stitch_coords(i).shift_in_x - min_shift_x + 1;

        % get coordinates for the lower right corner:
        scan.stitch_coords(i).corner_lr_row = scan.stitch_coords(i).corner_ul_row + image_height - 1;
        scan.stitch_coords(i).corner_lr_column = scan.stitch_coords(i).corner_ul_column + image_width - 1;

    end

end