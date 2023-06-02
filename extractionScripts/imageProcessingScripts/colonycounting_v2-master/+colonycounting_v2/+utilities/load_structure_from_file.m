function structure = load_structure_from_file(path_to_file)

    structure = load(path_to_file);
    
    % get fields:
    field_names = fieldnames(structure);

    % if there is only one field:
    if numel(field_names) == 1
        
        % load that field:
        structure = structure.(field_names{1});
    
    % otherwise:
    else
        
        % throw an error:
        error('The file contains more than 1 structure.');
        
    end
    
end