function image = add_boundary_to_image(image, boundary, color, line_width)

    % rearrange boundary to format preferred by insertShape:
    boundary_formatted = reshape(boundary', 1, []);

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
    
    % add to image:
    image = insertShape(image, 'Polygon', boundary_formatted, 'Color', color, 'LineWidth', line_width);

end