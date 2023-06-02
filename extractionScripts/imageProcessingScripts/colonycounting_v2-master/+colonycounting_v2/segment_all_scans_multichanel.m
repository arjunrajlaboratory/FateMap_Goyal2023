 function segment_all_scans_multichanel(varargin)

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
            %%% (or load it if it already exists). Also, we want to ask the
            %%% user how different wavelengths do they want to store
            %%% segments for
            
            % get name and path for file to store boundaries:
            file_name_boundaries = sprintf('Segment_Info_%s.mat', stitch_info.name_scan);
            
            % check for existence of boundaries file:
            if exist(fullfile(path_scan, file_name_boundaries), 'file') == 2
                
                % load the boundaries:
                boundaries = colonycounting_v2.utilities.load_structure_from_file(fullfile(path_scan, file_name_boundaries));
                
            % if the boundaries file does not already exist:
            else
                
                %Get the wavelegths for default groups
                default_groups = extractfield(stitch_info.images, 'wavelength');
                
                % ask user to assign a channel to each wavelength:
                prompt = {'How many seperate gropups would you like to segment?'};
                default = {'1'};
                
                nGroups = inputdlg(prompt, 'Enter number of groups', [1 200], default);
                nGroups = str2double(nGroups);
                %Ask user to assign names to the groups:
                prompt = sprintfc('Group%d', 1:nGroups);
                groups = inputdlg(prompt, 'Enter name of each segment group', [1 200], default_groups);
    
                % create structure to store boundaries:
                boundaries = struct;
                
                %First segment the well. 
                    boundaries.well.stitch_original.status = [];
                    boundaries.well.stitch_original.coordinates_boundary = [];
                    boundaries.well.stitch_original.coordinates_mask = [];
                
                    boundaries.well.stitch_small.status = [];
                    boundaries.well.stitch_small.coordinates_boundary = [];
                    boundaries.well.stitch_small.coordinates_mask = [];
                    
                    % segment the well:
                    instructions_well = 'well';
                    boundaries.well.stitch_small = colonycounting_v2.segment_all_scans.gui_to_segment_a_stitch(stitch_small, boundaries.well.stitch_small, instructions_well);

                    % save boundaries:
                    save(fullfile(path_scan, file_name_boundaries), 'boundaries');

                    %Next, for each group, add segments to boundaries
                
                for k=1:nGroups

                    boundaries.(groups{k}).colonies.stitch_original.status = [];
                    boundaries.(groups{k}).colonies.stitch_original.coordinates_boundary = [];
                    boundaries.(groups{k}).colonies.stitch_original.coordinates_mask = [];
                
                    boundaries.(groups{k}).colonies.stitch_small.status = [];
                    boundaries.(groups{k}).colonies.stitch_small.coordinates_boundary = [];
                    boundaries.(groups{k}).colonies.stitch_small.coordinates_mask = [];
                    
                    % segment the colonies:
                    instructions_colonies = sprintf('Segment colonies for group %d', k);
                    boundaries.(groups{k}).colonies.stitch_small = colonycounting_v2.segment_all_scans.gui_to_segment_a_stitch(stitch_small, boundaries.(groups{k}).colonies.stitch_small, instructions_colonies);

                    % save boundaries:
                    save(fullfile(path_scan, file_name_boundaries), 'boundaries');
                    

                    %%% Next, we want to determine the cells in each boundary. 
            
                    % get boundary coords in reference frame of original stitch:
                    boundaries = colonycounting_v2.segment_all_scans.scale_coords_boundary_up_multichanel(boundaries, groups{k}, stitch_info.scale_rows, stitch_info.scale_columns);
            
                    % determine cells in the well/colonies:
                    cells.well = colonycounting_v2.segment_all_scans.get_cells_within_boundaries(cells.all, boundaries.well);
                    cells.(groups{k}).colonies = colonycounting_v2.segment_all_scans.get_cells_within_boundaries(cells.all, boundaries.(groups{k}).colonies);


                    % save boundaries:
                    boundaries.cells.well = cells.well;
                    boundaries.cells.(groups{k}).colonies = cells.(groups{k}).colonies;
                    save(fullfile(path_scan, file_name_boundaries), 'boundaries');


                end
                
                %%% Next, we want to overlay the results on the small stitch
                %%% and save. Note that we overlay the results specifically on
                %%% the SMALL stitch because the necessary MATLAB functions do
                %%% not like very large images. 

                % overlay on stitch:
                stitch_small_annotated = colonycounting_v2.segment_all_scans.overlay_on_image_multichannel(stitch_small, cells, groups);

                % save the annotated images:
                imwrite(stitch_small_annotated, fullfile(path_scan, sprintf('Segment_Small_%s.tif', name_scan)));
            end
           
        end
        
    end

end