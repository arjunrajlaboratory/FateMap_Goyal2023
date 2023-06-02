function stitch_a_wavelength(scan, images)
       
    %%% First, we need to load all of the images:

    % load an image to get the image dimensions:
    temp = readmm(fullfile(scan.path_folder, images.list_images(1).name));
    image_height = temp.height;
    image_width = temp.width;

    % get number of positions:
    num_positions = numel(images.list_images);

    % create empty array to store the positions as a z-stack:
    image_stack = zeros(image_height, image_width, num_positions, 'uint16');

    % for each position:
    for j = 1:num_positions

        % get the name of the image:
        temp_image_name = sprintf('%s_%s_s%s_t1.TIF', scan.name_scan, images.wavelength, num2str(j));

        % load the image into the stack:
        image_stack(:,:,j) = im2uint16(imread(fullfile(scan.path_folder, temp_image_name)));

    end

    %%% Next, we need to add each position to the stitched image.

    % get size of stitched image:
    stitch_height = extractfield(scan, 'stitch_coords');
    stitch_height = max(extractfield(stitch_height{1}, 'corner_ul_row')) + image_height - 1;
    stitch_width = extractfield(scan, 'stitch_coords');
    stitch_width = max(extractfield(stitch_width{1}, 'corner_ul_column')) + image_width - 1;

    % create an empty array to store stitched image:
    image_stitched = zeros(stitch_height, stitch_width, 'uint16');

    % for each position (starting at the end and counting down):
    for j = num_positions:-1:1

        % get the image slice and convert to double:
        temp_image_slice = image_stack(:,:,j);

        % set 
        image_stitched(...
            scan.stitch_coords(j).corner_ul_row:scan.stitch_coords(j).corner_lr_row,...
            scan.stitch_coords(j).corner_ul_column:scan.stitch_coords(j).corner_lr_column ...
            ) = temp_image_slice;

    end

    % create a downsized version of the stitch:
    image_stitched_small = imresize(image_stitched, [image_height image_width]);

    % get the scale factor:
    scan.scale_rows = size(image_stitched, 1) / image_height;
    scan.scale_columns = size(image_stitched, 2) / image_width;

    % set file name:
    file_name = fullfile(scan.path_folder, sprintf('Stitch_Original_%s_%s.mat', scan.name_scan, images.channel));
    file_name_small = fullfile(scan.path_folder, sprintf('Stitch_Small_%s_%s.tif', scan.name_scan, images.channel));

    % set images to save:
    image_to_save = image_stitched;
    image_to_save_small = image_stitched_small;

    % save small image as tif:
    imwrite(image_to_save_small, file_name_small);

    % save large image as mat:
    save(file_name, 'image_to_save', '-v7.3');

    % save stitch info:
    save(fullfile(scan.path_folder, sprintf('Stitch_Info_%s.mat', scan.name_scan)), 'scan');

end