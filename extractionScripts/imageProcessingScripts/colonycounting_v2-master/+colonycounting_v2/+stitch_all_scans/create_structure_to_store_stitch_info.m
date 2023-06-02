function stitch_info = create_structure_to_store_stitch_info(paths)

    % get number of folders:
    num_folders = numel(paths);

    % create structure to store stitch info:
    [stitch_info(1:num_folders).path_folder] = deal('');
    [stitch_info(1:num_folders).scan_info] = deal(struct);

    % for each folder:
    for i = 1:num_folders

        % get path to folder:
        path_folder = paths{i};

        % get a list of all images in the folder:
        list_image_names = dir(fullfile(path_folder, 'Scan*.tif'));

        % get a list of all the scans in the folder:
        list_scan_names = extractfield(list_image_names, 'name');
        list_scan_names = unique(cellfun(@(x) x(1:7), list_scan_names, 'UniformOutput', false));

        % get number of scans:
        num_scans = numel(list_scan_names);
        
        % create structure to store scans:
        [scan_info(1:num_scans).path_folder] = deal('');
        [scan_info(1:num_scans).name_scan] = deal('');
        [scan_info(1:num_scans).images] = deal(struct);

        % for each scan:
        for j = 1:num_scans

            % get name of scan:
            name_scan = list_scan_names{j};

            % get list of images for the scan:
            list_image_names_scan = colonycounting_v2.utilities.get_structure_results_containing_string(list_image_names, 'name', name_scan);

            % get list of wavelengths for the scan:
            list_wavelengths = extractfield(list_image_names_scan, 'name');
            list_wavelengths = unique(cellfun(@(x) x(9:10), list_wavelengths, 'UniformOutput', false));
            
            % get number of wavelengths:
            num_wavelengths = numel(list_wavelengths);

            % create structure to store images:
            [images(1:num_wavelengths).wavelength] = deal('');
            [images(1:num_wavelengths).channel] = deal('');  
            [images(1:num_wavelengths).list_images] = deal('');
            
            % for each wavelength:
            for k = 1:num_wavelengths
                
                % get wavelength:
                wavelength = list_wavelengths{k};
                
                % get list of images for the wavelength:
                list_images = ...
                    colonycounting_v2.utilities.get_structure_results_containing_string(...
                    list_image_names_scan, 'name', list_wavelengths{k});
                
                % save:
                images(k).wavelength = wavelength;
                images(k).list_images = list_images;
                
            end
      
            % save:
            scan_info(j).path_folder = path_folder;
            scan_info(j).name_scan = name_scan;
            scan_info(j).images = images;
            
        end

        % save:
        stitch_info(i).path_folder = path_folder;
        stitch_info(i).scan_info = scan_info;
        
    end

end