% This script formats raw timelapse data off the Incucyte and converts it
% into 1 stack for each position. It was created specifically by Lauren Beck (for images
% from Organoids Round 170), so it probably cannot be generalized as is to any
% data from the Incucyte. 

% ask user for folder containing raw images:
path_images_raw = uigetdir;

% ask user where to put re-formatted images:
path_images_formatted = uigetdir;

% ask user for prefix used on each image:
prefix = inputdlg('What is the prefix used for each image? This is the part of the image names specific to the experiment.');

% get list of raw images:
list_raw_images = dir(fullfile(path_images_raw, sprintf('%s*.tif', prefix{1})));
list_raw_images = extractfield(list_raw_images, 'name');

% get list of positions:
list_positions = unique(cellfun(@(x) x(1:strfind(x, 'y')-6), list_raw_images, 'UniformOutput', false));

% for each position:
for i = 1:numel(list_positions)
   
    % get list of timepoints (raw images) for that position:
    list_raw_images_one_position = list_raw_images(contains(list_raw_images, list_positions{i}));
     
    % get the image name for saving:
    image_name_and_path = fullfile(path_images_formatted, [list_positions{i} '.tif']);
    
    % for all other timepoints:
    for j = 1:numel(list_raw_images_one_position)
        
        % load the image:
        image_slice = imread(fullfile(path_images_raw, list_raw_images_one_position{j}));
        
        % save the slice:
        if j == 1
            imwrite(image_slice, image_name_and_path);
        else
            imwrite(image_slice, image_name_and_path, 'writemode', 'append');
        end
        
    end
    
end