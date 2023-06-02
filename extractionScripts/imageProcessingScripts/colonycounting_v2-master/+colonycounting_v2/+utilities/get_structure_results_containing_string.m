function [structure_match, rows_output] = get_structure_results_containing_string(structure, field, string)

    % get the rows of whose values match the string:
    rows = ~cellfun(@isempty, strfind(extractfield(structure, field), string));

    % subset the structure:
    structure_match = structure(rows, 1);
    
    % return rows if prompted:
    if nargout > 1

        rows_output = rows;

    end

end