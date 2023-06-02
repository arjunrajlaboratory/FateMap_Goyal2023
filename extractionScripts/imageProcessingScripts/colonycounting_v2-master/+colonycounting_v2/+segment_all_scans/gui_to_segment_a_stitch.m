function boundaries = gui_to_segment_a_stitch(stitch, boundaries, instructions)

    % find the location of the sound bite:
    %path_notification = which('hummus.mp3');
    
    % load the sound bite to play:
    %[notification_sound, notification_sample_rate] = audioread(path_notification);
    
    % set default thresholds to use:
    threshold_min = 0.0;
    threshold_max = 1.0;

    % create GUI:
    handles = create_GUI;
    
    % make the figure larger;
    set(handles.figure, ...
        'Units', 'normalized', ...
            'Position', [0, 0, .5, .9]);
        
    % move gui to the center of the screen:
    movegui(handles.figure, 'center');
    
    % view image:
    view_image;
    
    % function to create GUI:
    function handles = create_GUI
        
        % set sizes:
        margin = 0.01;
        image_height = 0.6;
        
        figure_start_x = 0;
        figure_start_y = 0;
        figure_width = 1;
        figure_height = 1;
        
        width_full = figure_width - 2*margin;
        width_half = (width_full - margin)/2;
        height_all = (figure_height - image_height - 6*margin)/4;
        
        image_start_x = margin;
        image_start_y = margin;
        image_width = width_full;
        
        contrast_lower_start_x = margin;
        contrast_lower_start_y = image_start_y + image_height + margin;
        contrast_lower_width = width_full;
        contrast_lower_height = height_all;
        
        contrast_upper_start_x = margin;
        contrast_upper_start_y = contrast_lower_start_y + height_all + margin;
        contrast_upper_width = width_full;
        contrast_upper_height = height_all;
        
        add_start_x = margin;
        add_start_y = contrast_upper_start_y + height_all + margin;
        add_width = width_half;
        add_height = height_all;
        
        done_start_x = add_start_x + width_half + margin;
        done_start_y = add_start_y;
        done_width = width_half;
        done_height = height_all;
        
        instructions_start_x = margin;
        instructions_start_y = done_start_y + height_all + margin;
        instructions_width = width_full;
        instructions_height = height_all;
        
        % create figure:
        handles.figure = figure('Units', 'normalized', ...
            'Position', [figure_start_x, figure_start_y, figure_width, figure_height]);
        
        % add image:
        handles.image = axes('Units', 'normalized', ...
            'Position', [margin, margin, image_width, image_height]);
        
        % add button to adjust lower end of contrast:
        handles.contrast_lower = uicontrol('Style', 'slider', ...
            'Units', 'normalized', ...
            'Position', [contrast_lower_start_x, contrast_lower_start_y, contrast_lower_width, contrast_lower_height], ...
            'Callback', @callback_contrast_lower, ...
            'min', 0, 'max', 1, 'Value', 0);
        
        % add button to adjust upper end of contrast:
        handles.contrast_upper = uicontrol('Style', 'slider', ...
            'Units', 'normalized', ...
            'Position', [contrast_upper_start_x, contrast_upper_start_y, contrast_upper_width, contrast_upper_height], ...
            'Callback', @callback_contrast_upper, ...
            'min', 0, 'max', 1, 'Value', 1);
        
        % add button to add a segmentation:
        handles.add = uicontrol('Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [add_start_x, add_start_y, add_width, add_height], ...
            'String', 'Add', ...
            'FontSize', 20, ...
            'Callback', @callback_add);
        
        % add button to be done:
        handles.done = uicontrol('Style', 'pushbutton', ...
            'Units', 'normalized', ...
            'Position', [done_start_x, done_start_y, done_width, done_height], ...
            'String', 'Done', ...
            'FontSize', 20, ...
            'Callback', @callback_done);
        
        % add instructions:
        uicontrol('Style', 'text', ...
            'Units', 'normalized', ...
            'Position', [instructions_start_x, instructions_start_y, instructions_width, instructions_height], ...
            'FontSize', 26, ...
            'String', sprintf('Segment the %s. Add segmentations by clicking "Add", drawing on the image, and double-clicking within the annotation. Delete segmentations by right clicking on the segmentation to delete. Adjust the image contrast using the sliders.', instructions));
        
    end

    % function to view the image:
    function view_image
        
        % display the image:
        imshow(imadjust(stitch, [threshold_min, threshold_max]), 'Parent', handles.image);
        
        % plot boundaries:
        plot_boundaries;
        
        % establish callback for deleting segmentation:
        set(gcf, 'WindowButtonDownFcn', @callback_delete);
        
        % have program wait:
        uiwait(handles.figure);
        
    end

    % function to plot the boundaries:
    function plot_boundaries
       
        % turn on the hold:
        hold(handles.image, 'on'); 
        
        % for each boundary:
        for i = 1:numel(boundaries)
            
            % if boundary should be displayed:
            if strcmp(boundaries(i).status, 'keep')
                
                % get coordinates:
                temp_coords = boundaries(i).coordinates_boundary;
                
                % add first coordinate to end of coordinates (so boundary
                % plots closed):
                temp_coords(end+1, :) = temp_coords(1,:);
           
            	% plot coordinates on image:
                plot(temp_coords(:,1), temp_coords(:,2), 'red');
            
            end
            
        end
        
        % turn off the hold:
        hold(handles.image, 'off'); 
        
    end

    % function to adjust lower bound of contrast:
    function callback_contrast_lower(~,~)
        
        % get the slider value:
        threshold_min = get(handles.contrast_lower, 'Value');
        
        % show the image:
        view_image
        
    end

    % function to adjust upper bound of contrast:
    function callback_contrast_upper(~,~)
        
        % get the slider value:
        threshold_max = get(handles.contrast_upper, 'Value');
        
        % show the image:
        view_image;
        
    end

    % callback to add a segmentation:
    function callback_add(~, ~)
            
        % allow user to draw on the stitch:
        handle_boundary = imfreehand(handles.image, 'Closed', 'true');

        % create function that will constrain the user from drawing
        % boundaries outside the image:
        constrain_function = makeConstrainToRectFcn('imfreehand', [1 size(stitch,1)], [1 size(stitch, 2)]); 
        setPositionConstraintFcn(handle_boundary, constrain_function);

        % wait for user to double-click:
        wait(handle_boundary);

        % make a sound when the user has finished a boundary:
        %sound(notification_sound, notification_sample_rate);
        
        % get status of boundary:
        temp.status = 'keep';
        
        % get coordinates of boundary:
        temp.coordinates_boundary = getPosition(handle_boundary);

        % round boundary coordinates:
        temp.coordinates_boundary = round(temp.coordinates_boundary);

        % convert boundary coordinates to mask coordinates:
        mask = poly2mask(temp.coordinates_boundary(:,1), temp.coordinates_boundary(:,2), size(stitch, 1), size(stitch, 2));
        [temp.coordinates_mask(:,1), temp.coordinates_mask(:,2)] = find(mask == 1);

        % save coords:
        boundaries = colonycounting_v2.utilities.add_entry_to_structure(temp, boundaries);

        % view image:
        view_image;

    end

    % callback to delete a segmentation:
    function callback_delete(~, ~)
        
        % if the click was a right click:
        if strcmp(get(gcf, 'SelectionType'), 'alt')
           
            % ask user to click on image:
            [point] = get(gca, 'CurrentPoint');

            % round point:
            point = fliplr(round(point(1, 1:2)));
            
            % get coordiantes within radius of point:
            radius = 10;
            point_radius_x = point(1)-radius:point(1)+radius;
            point_radius_y = point(2)-radius:point(2)+radius;
            [point_radius_x, point_radius_y] = meshgrid(point_radius_x, point_radius_y);
            point_radius = cat(2, point_radius_x, point_radius_y);
            point_radius = reshape(point_radius, [], 2);

            % for each current boundary point:
            for i = 1:numel(boundaries)

                % if any point(s) falls within mask:
                if any(ismember(point_radius, boundaries(i).coordinates_mask, 'rows'))

                    % set that boundaries status to remove:
                    boundaries(i).status = 'remove';
                    
                    % make a sound when the user has finished a boundary:
                    %sound(notification_sound, notification_sample_rate);

                end

            end

            % view the image:
            view_image;
            
        end
        
    end

    % callback to be done:
    function callback_done(~,~)
        
        % remove any boundaries the user wanted to delete:
        boundaries = colonycounting_v2.utilities.get_structure_results_matching_string(boundaries, 'status', 'keep');
        
        % close the GUI:
        close(handles.figure);
        
    end

end