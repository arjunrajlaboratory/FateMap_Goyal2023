function image = add_cells_to_image(image, coordinates, color, marker_size)

    % format coordinates:
    coordinates_formatted = coordinates;

    % scale color to the image class:
    if isa(image, 'uint8')
        color = color * (2^8);
    elseif isa(image, 'uint16')
        color = color * (2^16);
    elseif isa(image, 'uint32')
        color = color * (2^32);
    elseif isa(image, 'uint64')
        color = color * (2^64);
    end

    % add coordinates to image:
    image = insertMarker(image, coordinates_formatted, 'circle', 'Color', color, 'Size', marker_size);

end