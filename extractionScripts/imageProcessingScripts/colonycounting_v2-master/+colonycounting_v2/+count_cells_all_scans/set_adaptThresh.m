function varargout = set_adaptThresh(stitch_small, stitch_original, scale_rows, scale_columns)
    %%% Need to update script to incorporate scaling option
    
    % set the default sensitivity and minimum nucleus size (in pixels):
    sensitivity = 0.3;
    minNucleusSize = 100;
    scaleOption = true;
    
    % set the default position:
    position_small.x = round(size(stitch_small, 2) / 2);
    position_small.y = round(size(stitch_small, 1) / 2);
    position_small.width = 100;
    position_small.height = 100;
    
    % create the figure:
    handles = create_GUI;

    % update the display:
    update_display;
    
    % make the figure visible:
    handles.figure.Visible = 'on';
    
    % have the function wait:
    uiwait(handles.figure);
    
    % output the sensitivity:
    varargout{1} = sensitivity;
    varargout{2} = minNucleusSize;
    varargout{3} = scaleOption;
    
    % function to create the GUI:
    function handles = create_GUI

        % set the margin size:
        margin = 0.001;

        % set the dimensions:
        width_full = 1 - (2 * margin);
        width_half = (1 - (3 * margin)) / 2;
        width_quarter = (1 - (5 * margin)) / 4;
        
        x_axes_full = margin;
        x_axes_zoomed = margin + width_half + margin;
        x_input_full = margin;
%         x_input_scale = margin + width_quarter + margin;
        x_input_zoomed = margin + width_half + margin;
        x_input_zoomed2 = margin + width_half + margin + width_quarter + margin;
        x_label_full = margin;
%         x_label_scale = margin + width_quarter + margin;
        x_label_zoomed = margin + width_half + margin;
        x_label_zoomed2 = margin + width_half + margin + width_quarter + margin;
        x_instructions = margin;

        height_instructions = 0.1;
        height_label = 0.025;
        height_input = 0.025;
        height_axes = 1 - height_instructions - height_label - height_input - (5 * margin);

        y_axes_full = margin;
        y_axes_zoomed = margin;
        y_input_full = margin + height_axes + margin;
        y_input_zoomed = margin + height_axes + margin;
        y_label_full = margin + height_axes + margin + height_input + margin;
        y_label_zoomed = margin + height_axes + margin + height_input + margin;
        y_instructions = margin + height_axes + margin + height_input + margin + height_label + margin;

        % create the figure window:
        handles.figure = figure('Visible', 'off', 'Units', 'normalized', 'Position', [0, 0, 0.8, 0.8]);

        % create axes for images:
        handles.axes_full = axes(handles.figure, 'Units', 'normalized', 'Position', [x_axes_full, y_axes_full, width_half, height_axes]);
        handles.axes_zoomed = axes(handles.figure, 'Units', 'normalized', 'Position', [x_axes_zoomed, y_axes_zoomed, width_half, height_axes]);
        
        % establish callback for clicking on the image:
        set(gcf, 'WindowButtonDownFcn', @callback_change_position);
        
        % create the input for box size:
        handles.input_full = uicontrol(handles.figure, 'Style', 'edit', 'Units', 'normalized', 'Position', [x_input_full, y_input_full, width_half, height_input], 'FontSize', 16);
        
        % establish callback for changing the box size:
        handles.input_full.Callback = @callback_change_box_size;
        
        % create the input for filter size:
        handles.input_zoomed = uicontrol(handles.figure, 'Style', 'edit', 'Units', 'normalized', 'Position', [x_input_zoomed, y_input_zoomed, width_quarter, height_input], 'FontSize', 16);
  
        % establish callback for changing the filter size:
        handles.input_zoomed.Callback = @callback_change_filter_size;
        
        % create the input for nucleus size:
        handles.input_zoomed2 = uicontrol(handles.figure, 'Style', 'edit', 'Units', 'normalized', 'Position', [x_input_zoomed2, y_input_zoomed, width_quarter, height_input], 'FontSize', 16);
  
        % establish callback for changing the minimum nucleus size:
        handles.input_zoomed2.Callback = @callback_change_minNucleusSize;
        
%         % create the input for scale option:
%         handles.input_scale = uicontrol(handles.figure, 'Style', 'togglebutton', 'Units', 'normalized', 'Position', [x_input_scale, y_input_zoomed, width_quarter, height_input], 'FontSize', 16);
%   
%         % establish callback for changing the minimum nucleus size:
%         handles.input_scale.Callback = @callback_change_scale;
%         
        % create the label for box size:
        handles.label_full = uicontrol(handles.figure, 'Style', 'text', 'Units', 'normalized', 'Position', [x_label_full, y_label_full, width_half, height_label], 'FontSize', 16);
        handles.label_full.String = 'size of field of view (units are meaningless)';

        % create the label for filter size:
        handles.label_zoomed = uicontrol(handles.figure, 'Style', 'text', 'Units', 'normalized', 'Position', [x_label_zoomed, y_label_zoomed, width_quarter, height_label], 'FontSize', 16);
        handles.label_zoomed.String = 'size of the gaussian filter';
        
        % create the label for minNucleusSize:
        handles.label_zoomed2 = uicontrol(handles.figure, 'Style', 'text', 'Units', 'normalized', 'Position', [x_label_zoomed2, y_label_zoomed, width_quarter, height_label], 'FontSize', 16);
        handles.label_zoomed2.String = 'minimum size of nucleus';

%         % create the label for the scale option:
%         handles.label_scale = uicontrol(handles.figure, 'Style', 'text', 'Units', 'normalized', 'Position', [x_label_scale, y_label_zoomed, width_quarter, height_label], 'FontSize', 16);
%         handles.label_scale.String = 'Option to scale';

        % create the instructions:
        handles.instructions = uicontrol(handles.figure, 'Style', 'text', 'Units', 'normalized', 'Position', [x_instructions, y_instructions, width_full, height_instructions], 'FontSize', 16);
        handles.instructions.String = {'1. Set the sensitivty of the adaptive threshold (must be 0-1) to use for counting cells. Also set a minimun nucleus size (in pixels) for exlcuding debris. Hit enter to update the display (it may take a few moments).', ...
            '2. Check that the values work well for a few different fields of view. You can change the field of view by clicking anywhere on the left image. You can also change the size of the field of the view (too large may slow down the code).', ...
            '3. Close the window when you are satisfied with the results for a particular gaussian size.'};
        
        % move gui to the center of the screen:
        movegui(handles.figure, 'center');

    end

    % function to update the display:
    function update_display

        % fill out the boxes:
        set(handles.input_zoomed, 'String', num2str(sensitivity));
        set(handles.input_zoomed2, 'String', num2str(minNucleusSize));
        set(handles.input_full, 'String', num2str(position_small.width));
        
        % convert the position coords from the small to original stitch size:
        position_original.width = position_small.width * scale_columns;
        position_original.height = position_small.height * scale_rows;
        position_original.x = position_small.x * scale_rows;
        position_original.y = position_small.y * scale_columns;

        % crop the original stitch:
        stitch_original_cropped = imcrop(stitch_original, [position_original.x, position_original.y, position_original.width, position_original.height]);

        % count cells in the cropped stitch:
        cells = colonycounting_v2.count_cells_all_scans.count_cells_in_scan_adaptThresh(stitch_original_cropped, sensitivity, minNucleusSize, scaleOption);

        % display small stitch:
        imshow(scale(stitch_small), 'Parent', handles.axes_full);

        % display the cropped original stitch:
        imshow(scale(stitch_original_cropped), 'Parent', handles.axes_zoomed);

        % plot the selected position:
        hold(handles.axes_full, 'on');
        rectangle(handles.axes_full, 'Position', [position_small.x, position_small.y, position_small.width, position_small.height], 'EdgeColor', 'red', 'LineWidth', 2);
        hold(handles.axes_full, 'off');

        % plot the cells:
        hold(handles.axes_zoomed, 'on');
        scatter(handles.axes_zoomed, cells(:,1), cells(:,2), 6, 'r', 'filled');
        hold(handles.axes_zoomed, 'off');
        
        % plot a red outline around the image:
        border = [1, 1; size(stitch_original_cropped, 1) - 1, 1; size(stitch_original_cropped, 1) - 1, size(stitch_original_cropped, 2) - 1; 1, size(stitch_original_cropped, 2) - 1; 1, 1];
        hold(handles.axes_zoomed, 'on');
        plot(handles.axes_zoomed, border(:,1), border(:,2), '-r', 'LineWidth', 2);
        hold(handles.axes_zoomed, 'off');

    end

    % function to change the position:
    function callback_change_position(~, ~)

        % get coordinates of selection:
        point = get(gca, 'CurrentPoint');

        % format coordinates:
        point = round(point(1,:));
        
        % if the coordinates are within the bounds of the full axes:
        if (point(1) > 1 && point(1) < size(stitch_small, 2)) && (point(2) > 1 && point(2) < size(stitch_small, 1))

            % update the position:
            position_small.x = point(1) - (position_small.width / 2);
            position_small.y = point(2) - (position_small.height / 2);

            % update the display:
            update_display;
        
        end

    end

    % function to change the filter size:
    function callback_change_filter_size(~, ~)
        
        % get the new sensitivity:
        sensitivity = str2double(handles.input_zoomed.String);
        
        % update the display:
        update_display;
        
    end

    % function to change the minimum nucleus size:
    function callback_change_minNucleusSize(~, ~)
        
        % get the new minNucleusSize:
        minNucleusSize = str2double(handles.input_zoomed2.String);
        
        % update the display:
        update_display;
        
    end

%     % function to toggle scale option:
%     function callback_change_scale(~, ~)
%         
%         % get toggle state
%         button_state = get(handles.input_scale,'Value');
% 
%         if button_state == get(handles.input_scale,'Max')
%             scaleOption = true;
%             set(handles.input_scale, 'String', 'true');
%         elseif scaleOption == get(handles.input_scale,'Min')
%             scaleOption = false;
%             set(handles.input_scale, 'String', 'false');
%         end
%         
%         % update the display:
%         update_display;
%         
%     end

% function to change the box size:
    function callback_change_box_size(~, ~)
        
        position_small.width = str2double(handles.input_full.String);
        position_small.height = str2double(handles.input_full.String);

        % update the display:
        update_display;
        
    end

end