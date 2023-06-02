function [structure_match, rows_output] = get_structure_results_matching_string(structure, field, string)

    % get the rows of whose values match the string:
    rows = find(strcmp(extractfield(structure, field), string) == 1);

    % subset the structure:
    structure_match = structure(1, rows);
    
    % return rows if prompted:
    if nargout > 1

        rows_output = rows;

    end

end