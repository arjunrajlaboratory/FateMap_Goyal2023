function stitch_info = get_scan_dimensions(stitch_info)

    % ask user if any scans have the same dimensions:
    question = 'Do any scans have the same sizes?';
    title = 'Set Scan Size';
    answer = questdlg(question, title, 'Yes, all do.', 'Yes, all in a folder do.', 'No, none do.', 'Yes, all do.');
     
    % if all scans have the same sizes:
    if strcmp(answer, 'Yes, all do.')
       
        % get scan size:
        [num_rows, num_columns] = colonycounting_v2.stitch_all_scans.get_scan_dimensions.enter_scan_dimensions('ALL FOLDERS', 'ALL SCANS');
        
    end
    
    % for each folder:
    for i = 1:numel(stitch_info)
       
        % if all scans in the folder have the same sizes:
        if strcmp(answer, 'Yes, all in a folder do.')

            % get scan size:
            [num_rows, num_columns] = colonycounting_v2.stitch_all_scans.get_scan_dimensions.enter_scan_dimensions(stitch_info(i).path_folder, 'ALL SCANS');

        end
        
        % for each scan:
        for j = 1:numel(stitch_info(i).scan_info)
           
            % if no scans have the same sizes:
            if strcmp(answer, 'No, none do.')

                % get channel names:
                [num_rows, num_columns] = colonycounting_v2.stitch_all_scans.get_scan_dimensions.enter_scan_dimensions(stitch_info(i).path_folder, stitch_info(i).scan_info(j).name_scan);

            end
            
            % save:
            stitch_info(i).scan_info(j).num_tiles_row = num_rows;
            stitch_info(i).scan_info(j).num_tiles_column = num_columns;
            stitch_info(i).scan_info(j).num_tiles = num_rows * num_columns;
            
            % for each wavelength:
            for k = 1:numel(stitch_info(i).scan_info(j).images)
                
                % get list of images for the wavelength:
                list_images = stitch_info(i).scan_info(j).images(k).list_images;
                
                % remove the Metamorph stitched image from the list of images:
                [~, rows_to_remove] = colonycounting_v2.utilities.get_structure_results_containing_string(list_images, 'name', num2str((num_rows*num_columns) + 1));
                list_images(rows_to_remove) = [];

                % save to stitch info structure:
                stitch_info(i).scan_info(j).images(k).list_images = list_images;
                
            end
            
        end
        
    end

end