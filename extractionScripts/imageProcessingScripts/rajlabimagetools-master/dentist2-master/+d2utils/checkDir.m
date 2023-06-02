function TF = checkDir(x)
    TF = false;
    if ~ischar(x)
        error('Input must be character not %s', class(x))
    elseif isempty(x)
        TF = true;
    elseif ~isfolder(x)
        error('Unable to find directory %s', x)
    else
        TF = true;
    end
end