function [shift_column, shift_row] = get_shift_overlap(folder, name_scan)

    % ask user for overlap:
    prompts = {folder, 'Column X Overlap', 'Column Y Overlap', 'Row X Overlap', ' Row Y Overlap'};
    defaults = cell(numel(prompts), 1);
    defaults(:,1) = {''};
    defaults{1} = name_scan;
    overlap = inputdlg(prompts, 'Enter Overlap', [1 200], defaults);
    
    % get overlaps:
    overlap_column_x = str2double(overlap{2});
    overlap_column_y = str2double(overlap{3});
    overlap_row_x = str2double(overlap{4});
    overlap_row_y = str2double(overlap{5});
    
    % divide all overlaps by 2 (since user enters the amount of
    % overlap, and the amount you need to shift the pixels by is HALF of
    % that distance):
    overlap_column_x = overlap_column_x / 2;
    overlap_column_y = overlap_column_y / 2;
    overlap_row_x = overlap_row_x / 2;
    overlap_row_y = overlap_row_y / 2;
    
    % take max between overlap and 0 (in case overlap was entered as 0 and
    % 0/2 is infinity):
    overlap_column_x = max(overlap_column_x, 0);
    overlap_column_y = max(overlap_column_y, 0);
    overlap_row_x = max(overlap_row_x, 0);
    overlap_row_y = max(overlap_row_y, 0);
    
    % split overlap into format used by cpselect: 
    shift_column.column = overlap_column_x;
    shift_column.row = overlap_column_y;
    shift_row.column = overlap_row_x;
    shift_row.row = overlap_row_y;
    
end