function image_density = create_density_image(image_original, window)

    % create empty array:
    image_density = zeros(size(image_original));
    
    % get image dimensions:
    image_height = size(image_original, 1);
    image_width = size(image_original, 2);

    % for each row:
    for i = 1:image_height

        % for each column:
        for j = 1:image_width

            % get x and y coordinates of the window to use:
            y_min = max(i - window, 1);
            y_max = min(i + window, image_height);
            x_min = max(j - window, 1);
            x_max = min(j + window, image_width);

            % get the number of non-zero elements (cells) within the window:
%             image_density(i,j) = nnz(image_original(y_min:y_max, x_min:x_max));
            image_density(i,j) = sum(sum(image_original(y_min:y_max, x_min:x_max)));
            
        end

    end

end