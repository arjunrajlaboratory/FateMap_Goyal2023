 function segment_all_scans(varargin)

    %%% First, we need to get all the path(s) to the raw images. 

    paths = colonycounting_v2.utilities.get_paths_to_data(varargin);
    
    %%% Next, we need to segment the wells.
    
    % for each path:
    for i = 1:numel(paths)
       
        % get path to scan:
        path_scan = paths{i};
        
        % get list of stitch infos (scans) in the folder:
        list_stitch_info = dir(fullfile(path_scan, 'Stitch_Info_*.mat'));
        
        % for each stitch info (scan) in the folder:
        for j = 1:numel(list_stitch_info)
           
            %%% First, we want to load the stitch and stitch info.
            
            % load the stitch info:
            stitch_info = colonycounting_v2.utilities.load_structure_from_file(fullfile(path_scan, list_stitch_info(j).name));
            
            % get the name of the scan:
            name_scan = stitch_info.name_scan;
            
            % load the cells:
            cells = colonycounting_v2.utilities.load_structure_from_file(fullfile(path_scan, sprintf('Cell_Info_%s.mat', name_scan)));
            
            % get the name of the small DAPI stitch:
            name_stitch_small = sprintf('Stitch_Small_%s_%s.tif', name_scan, 'dapi');
            
            % load the small stitch:
            stitch_small = readmm(fullfile(path_scan, name_stitch_small));
            stitch_small = stitch_small.imagedata;
            
            %%% Next, we want to create a structure to store the boundaries
            %%% (or load it if it already exists). Also, we want to guess
            %%% what the colony segmentations should be. 
            
            % get name and path for file to store boundaries:
            file_name_boundaries = sprintf('Segment_Info_%s.mat', stitch_info.name_scan);
            
            % check for existence of boundaries file:
            if exist(fullfile(path_scan, file_name_boundaries), 'file') == 2
                
                % load the boundaries:
                boundaries = colonycounting_v2.utilities.load_structure_from_file(fullfile(path_scan, file_name_boundaries));
                
            % if the boundaries file does not already exist:
            else
                
                % create structure to store boundaries:
                boundaries = struct;
                
                boundaries.well.stitch_original.status = [];
                boundaries.well.stitch_original.coordinates_boundary = [];
                boundaries.well.stitch_original.coordinates_mask = [];
                
                boundaries.well.stitch_small.status = [];
                boundaries.well.stitch_small.coordinates_boundary = [];
                boundaries.well.stitch_small.coordinates_mask = [];
                
                boundaries.colonies.stitch_original.status = [];
                boundaries.colonies.stitch_original.coordinates_boundary = [];
                boundaries.colonies.stitch_original.coordinates_mask = [];
                
                boundaries.colonies.stitch_small.status = [];
                boundaries.colonies.stitch_small.coordinates_boundary = [];
                boundaries.colonies.stitch_small.coordinates_mask = [];
                
            end
            
            %%% Next, we want the user to complete/review the
            %%% segmentations. Note that I save the boundaries periodically
            %%% to minimize the event that errors in the GUI will cause the
            %%% function to crash (and the user would lose their work).
            
            % segment the well:
            instructions_well = 'well';
            boundaries.well.stitch_small = colonycounting_v2.segment_all_scans.gui_to_segment_a_stitch(stitch_small, boundaries.well.stitch_small, instructions_well);
            
            % save boundaries:
            save(fullfile(path_scan, file_name_boundaries), 'boundaries');
            
            % ask the user if they want the computer to guess the colonies:
            answer = questdlg('Do you want the computer to guess the colony boundaries?', 'Colony Segmentation', 'Yes', 'No', 'No');
            
            % if the user wants to have the computer guess the colony
            % boundaries:
            if strcmp(answer, 'Yes')
            
                % guess the colony boundaries: 
                boundaries.colonies.stitch_small = colonycounting_v2.segment_all_scans.guess_colonies(cells.all.stitch_small, boundaries.colonies.stitch_original, stitch_small, boundaries.well.stitch_small);
            
            end
            
            % segment the colonies:
            instructions_colonies = 'colonies';
            boundaries.colonies.stitch_small = colonycounting_v2.segment_all_scans.gui_to_segment_a_stitch(stitch_small, boundaries.colonies.stitch_small, instructions_colonies);

            % save boundaries:
            save(fullfile(path_scan, file_name_boundaries), 'boundaries');
            
            %%% Next, we want to determine the cells in each boundary. 
            
            % get boundary coords in reference frame of original stitch:
            boundaries = colonycounting_v2.segment_all_scans.scale_coords_boundary_up(boundaries, stitch_info.scale_rows, stitch_info.scale_columns);
            
            % determine cells in the well/colonies:
            cells.well = colonycounting_v2.segment_all_scans.get_cells_within_boundaries(cells.all, boundaries.well);
            cells.colonies = colonycounting_v2.segment_all_scans.get_cells_within_boundaries(cells.all, boundaries.colonies);

            %%% Next, we want to overlay the results on the small stitch
            %%% and save. Note that we overlay the results specifically on
            %%% the SMALL stitch because the necessary MATLAB functions do
            %%% not like very large images. 
                       
            % overlay on stitch:
            stitch_small_annotated = colonycounting_v2.segment_all_scans.overlay_on_image(stitch_small, cells);
            
            % save boundaries:
            boundaries.cells.well = cells.well;
            boundaries.cells.colonies = cells.colonies;
            save(fullfile(path_scan, file_name_boundaries), 'boundaries');
            
            % save the annotated images:
            imwrite(stitch_small_annotated, fullfile(path_scan, sprintf('Segment_Small_%s.tif', name_scan)));

        end
        
    end

end