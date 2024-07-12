function plot_point_pages(max_points, points)
    figure;
    num_pages = max_points;
    for k = 1:num_pages
        subplot(ceil(sqrt(num_pages)), ceil(sqrt(num_pages)), k);
        hold on;
        imagesc(points(:, :, k));
        colorbar;
        clim([0, 1e7]);
        xlabel('Pixel Index (Column)');
        ylabel(['Value of Point ' num2str(k)]);
        title(['Points in Bucket ' num2str(k)]);
        axis equal;
        axis tight;
        hold off;
    end
end

% Step 1: Initialize points for each pixel
% Example initialization for illustration purposes
image_size = [180, 180]; % Replace with your image size
max_points = 6; % Ordinarily the maximum length would have been computed during construction but good enough
points = zeros(image_size(1), image_size(2), max_points);
for i = 1:image_size(1)
    for j = 1:image_size(2)
        num_points = randi([2, 6]); % Random number of points between 2 and 6
        dummy_peaks = sort(rand(1, num_points)) .* 1e7; % Random values for the points
        points(i, j, 1:num_points) = dummy_peaks; % Embed the actual peaks into a nicer version
    end
end

plot_point_pages(max_points, points)

% Difference metric
difference = @(p1, p2) abs(p1 - p2);

% Initialize ordering
ordering = cell(image_size);
for i = 1:image_size(1)
    for j = 1:image_size(2)
        num_points = length(points(i, j, :));
        ordering{i, j} = 1:num_points; % Initial ordering is the original order
    end
end

% Apply the best ordering to the points
for i = 1:image_size(1)
    for j = 1:image_size(2)
        points(i, j, :) = points(i, j, best_ordering{i, j});
    end
end

% Step 4: Visualize the points
% Step 4: Visualize the points in paged subplots
plot_point_pages(max_points, points)