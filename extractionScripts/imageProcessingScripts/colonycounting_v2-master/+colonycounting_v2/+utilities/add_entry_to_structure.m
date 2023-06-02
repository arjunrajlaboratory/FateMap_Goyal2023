function structure = add_entry_to_structure(entry, structure)

    % get field names:
    field_names = fieldnames(structure);

    % if the first field of the structure is empty:
    if isempty(extractfield(structure, field_names{1}))
        
        % use the entry to define the structure:
        structure = entry;
        
    % otherwise:
    else
        
        % add entry to end of structure:
        structure = [structure, entry];
        
    end

end