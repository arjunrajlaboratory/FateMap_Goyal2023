function TF = checkFile(x)
    TF = false;
    if ~ischar(x)
        error('Input must be character not %s', class(x))
    elseif isempty(x)
        TF = true;
    elseif ~isfile(x)
        error('Unable to find file %s', x)
    else
        TF = true;
    end
end