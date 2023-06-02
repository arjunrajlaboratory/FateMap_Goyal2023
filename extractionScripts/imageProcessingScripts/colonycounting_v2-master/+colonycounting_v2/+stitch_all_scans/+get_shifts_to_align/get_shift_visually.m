function shift = get_shift_visually(image_left, image_right)

    % register the image using cp select:
    [moving, fixed] = cpselect(scale(image_left), scale(image_right), 'Wait', true);
       
    % get the distance between each moving and fixed coords:
    distance_to_transform = moving - fixed;

    % get the median distance (for when the user set multiple points):
    distance_to_transform = median(distance_to_transform, 1);

    % round the distance to a whole number (since you cannot have half a pixel):
    distance_to_transform = round(distance_to_transform);
    
    % save:
    shift.row = distance_to_transform(1);
    shift.column = distance_to_transform(2);

end