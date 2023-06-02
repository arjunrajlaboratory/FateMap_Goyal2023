function image_annotated = overlay_on_image(image, cells)

    % create rgb-version of image:
    image_annotated = repmat(image, 1, 1, 3);

    % set color:
    color = [0 0 1];
    
    % add cells to image:
    image_annotated = colonycounting_v2.utilities.add_cells_to_image(image_annotated, cells, color, 1);

end