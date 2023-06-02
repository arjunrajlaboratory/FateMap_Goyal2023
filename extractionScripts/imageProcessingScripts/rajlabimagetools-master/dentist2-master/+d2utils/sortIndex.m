function idx = sortIndex(arrayIn, direction)
    [~, idx] = sort(arrayIn, direction, 'MissingPlacement', 'last');
end