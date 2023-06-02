function [colony_centroids_x, colony_centroids_y] = get_colony_centroids(image)

    % blur the density image:
    image_blur = imgaussfilt(image, 10);

    % get regionl max:
    regional_max = imregionalmax(image_blur);

    % get centroids of regional max (centroids of colonies):
    colony_centroids = regionprops(regional_max, 'Centroid');

    % get number of colonies:
    num_colonies = numel(colony_centroids);

    % format centroids:
    colony_centroids_x = zeros(num_colonies, 1);
    colony_centroids_y = zeros(num_colonies, 1);
    for i = 1:num_colonies
       colony_centroids_x(i) = round(colony_centroids(i).Centroid(1)); 
       colony_centroids_y(i) = round(colony_centroids(i).Centroid(2));
    end

end