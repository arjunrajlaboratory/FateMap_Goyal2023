function paths = get_paths_to_data(varargin)

    % if no paths are input:
    if isempty(varargin{1})
        
        % set the current working directory as the path to the data:
        paths = {pwd};
        
    % otherwise:
    else
        
        % use the supplied cell as the path(s) to the data:
        paths = varargin{1}{1};
        
    end

end