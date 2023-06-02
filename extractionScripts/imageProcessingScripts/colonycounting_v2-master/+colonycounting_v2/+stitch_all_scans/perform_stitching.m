function perform_stitching(stitch_info)

    % for each folder:
    for i = 1:numel(stitch_info)
       
        % for each scan:
        for j = 1:numel(stitch_info(i).scan_info)

            % display status:
            colonycounting_v2.utilities.display_status('Stitching', stitch_info(i).scan_info(j).name_scan, stitch_info(i).path_folder);

            % for each wavelength:
            for k = 1:numel(stitch_info(i).scan_info(j).images)

                % stitch the wavelength:
                colonycounting_v2.stitch_all_scans.perform_stitching.stitch_a_wavelength(stitch_info(i).scan_info(j), stitch_info(i).scan_info(j).images(k));
            
            end

        end
        
    end
    
end