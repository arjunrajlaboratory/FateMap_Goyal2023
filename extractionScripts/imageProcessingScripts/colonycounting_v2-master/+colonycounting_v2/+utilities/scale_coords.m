function coords = scale_coords(coords, scale_rows, scale_columns)

    % scale the row coordinates:
    coords(:,1) = round(coords(:,1) * scale_columns);
    
    % scale the column coordinates:
    coords(:,2) = round(coords(:,2) * scale_rows);

end

