function structure_sorted = sort_structure_based_on_field(structure, field_to_sort_by)

    % get field names:
    field_names = fieldnames(structure);
    
    % get field number to sort by:
    field_index = find(contains(field_names, field_to_sort_by));
    
    % convert to cell:
    cell = struct2cell(structure);
    
    % get dimensions:
    sz = size(cell);
    
    % resize:
    cell = reshape(cell, sz(1), []);

    % make each field a column:
    cell = cell';          

    % sort by first field "name":
    cell_sorted = sortrows(cell, field_index);

    % put back into original cell array format:
    cell_sorted = reshape(cell_sorted', sz);

    % convert to back to structure:
    structure_sorted = cell2struct(cell_sorted, field_names, 1);

end