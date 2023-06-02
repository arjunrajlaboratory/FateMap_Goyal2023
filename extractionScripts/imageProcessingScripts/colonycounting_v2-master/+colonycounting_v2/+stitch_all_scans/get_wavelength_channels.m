function stitch_info = get_wavelength_channels(stitch_info)

    % ask user if any scans have the same wavelengths:
    question = 'Do any scans have the same channel-wavelength assignments?';
    title = 'Channel and Wavelength Assignment';
    answer = questdlg(question, title, 'Yes, all do.', 'Yes, all in a folder do.', 'No, none do.', 'Yes, all do.');
     
    % if all scans have the same wavelengths:
    if strcmp(answer, 'Yes, all do.')
       
        % get channel names:
        channels = colonycounting_v2.stitch_all_scans.get_wavelength_channels.get_channels_for_each_wavelength('ALL FOLDERS', 'ALL SCANS', stitch_info(1).scan_info(1).images);
        
    end
    
    % for each folder:
    for i = 1:numel(stitch_info)
       
        % if all scans in the folder have the same wavelengths:
        if strcmp(answer, 'Yes, all in a folder do.')

            % get channel names:
            channels = colonycounting_v2.stitch_all_scans.get_wavelength_channels.get_channels_for_each_wavelength(stitch_info(i).path_folder, 'ALL SCANS', stitch_info(i).scan_info(1).images);

        end
        
        % for each scan:
        for j = 1:numel(stitch_info(i).scan_info)
           
            % if no scans have the same wavelengths:
            if strcmp(answer, 'No, none do.')

                % get channel names:
                channels = colonycounting_v2.stitch_all_scans.get_wavelength_channels.get_channels_for_each_wavelength(stitch_info(i).path_folder, stitch_info(i).scan_info(j).name_scan, stitch_info(i).scan_info(j).images);

            end
            
            % for each wavelength:
            for k = 1:numel(stitch_info(i).scan_info(j).images)
                
                % save:
                stitch_info(i).scan_info(j).images(k).channel = channels{k};
                
            end
            
        end
        
    end

end

