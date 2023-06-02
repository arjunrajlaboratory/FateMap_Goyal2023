function downsize_all_scans(varargin)

    %%% First, we need to get all the path(s) to the raw images. 
    
    paths = colonycounting_v2.utilities.get_paths_to_data(varargin);
    
    %%% Next, we want to ask the user by how much to downsize the images. 
    
    % get user input:
    downsize_factor = inputdlg('How much do you want to downsize the images?', 'Downsizing the Images');
    
    % convert the downsize factor to a number:
    downsize_factor = str2double(downsize_factor{1});
    
    % round the downsize factor to a whole number:
    downsize_factor = round(downsize_factor);
    
    %%% Next, we want to create downsized versions of the images. 
    
    % for each folder:
    for i = 1:numel(paths)
        
        % get name of folder with original images:
        folder_original = paths{i};
        
        % get name of folder to store downsized images:
        folder_downsized = fullfile(folder_original, 'downsized_images');
        mkdir(folder_downsized);
       
        % get list of images:
        list_images = dir(fullfile(folder_original, 'Scan*.tif'));
        
        % for each image:
        for j = 1:numel(list_images)
            
            % get the image name:
            image_name = list_images(j).name;
            
            % load the image:
            image = imread(fullfile(folder_original, image_name));
            
            % downsize the image:
            image_downsize = imresize(image, (1/downsize_factor));
            
            % save the image:
            imwrite(image_downsize, fullfile(folder_downsized, image_name));
            
        end
        
    end

end