function plot_point_pages(max_points, points)
    figure;
    num_pages = max_points;
    intervals = zeros(max_points, 3);
    for k = 1:num_pages
        data = squeeze(points(:, :, k));
        data(data == 0) = NaN;
        intervals(k, :) = [min(data(:)), max(data(:)), range(data(:))];
        subplot(ceil(sqrt(num_pages)), ceil(sqrt(num_pages)), k);
        imAlpha=ones(size(data));
        imAlpha(isnan(data))=0;
        set(subplot(ceil(sqrt(num_pages)), ceil(sqrt(num_pages)), k), 'Color', 'Red');
        hold on;
        imagesc(data, 'AlphaData', imAlpha);
        colorbar;
        xlabel('Pixel Index (Column)');
        ylabel(['Value of Point ' num2str(k)]);
        title(['Points in Bucket ' num2str(k)]);
        axis equal;
        axis tight;
        hold off;
    end

    disp(intervals)
    % for i=1:size(intervals, 1)
    %     for j=i+1:size(intervals, 1)
    %         if (intervals(i, 1) < intervals(j, 1) && intervals(i, 2) > intervals(j, 1)) || ...
    %            (intervals(i, 1) < intervals(j, 2) && intervals(i, 2) > intervals(j, 2)) || ...
    %            (intervals(i, 1) > intervals(j, 1) && intervals(i, 2) < intervals(j, 2)) || ...
    %            (intervals(i, 1) < intervals(j, 1) && intervals(i, 2) > intervals(j, 2))
    %             disp("Failure")
    %             disp(intervals(i, :))
    %             disp(intervals(j, :))
    %         end
    %     end
    % end

    figure;
    subplot(1, 2, 1)
    hold on;
    % Horizontal line cut
    startRow = 90;
    data = squeeze(points(startRow, :, :));
    data(data == 0) = NaN;
    colors = ['r', 'g', 'b', 'c', 'm', 'k'];
    for i=1:size(points, 3)
        color = colors(mod(i - 1, numel(colors)) + 1);
        scatter(data(:, i), 1:size(points, 1), 10, color, 'filled', 'MarkerEdgeColor', color)
    end
    set(gca, 'YDir','reverse')
    title("Horizontal")
    hold off;
    subplot(1, 2, 2)
    hold on;
    % Vertical line cut
    startCol = 90;
    data = squeeze(points(:, startCol, :));
    data(data == 0) = NaN;
    colors = ['r', 'g', 'b', 'c', 'm', 'k'];
    for i=1:size(points, 3)
        color = colors(mod(i - 1, numel(colors)) + 1);
        scatter(data(:, i), 1:size(points, 1), 10, color, 'filled', 'MarkerEdgeColor', color)
    end
    set(gca, 'YDir','reverse')
    title("Vertical")
    hold off;
end

f = waitbar(0, "Initializing");
max_points = 12;
points = locs_array;
image_size = size(points);

waitbar(1/4, f, "Removing Zeroes");
i = 1;
for k=1:max_points
    if all(locs_array(:, :, k) == 0)
        points(:, :, i) = [];
    else    
        i = i + 1;
        max_points = max_points - 1;
    end
end

% Flatten the points matrix into a single list
waitbar(2/4, f, "Flattening");
total_points = nnz(points(:) ~= 0);
all_points = zeros(total_points, 1);
pixel_indices = zeros(total_points, 2);

% Initialize counters
point_counter = 1;
% Loop through the image and gather points and their indices
for i = 1:image_size(1)
    for j = 1:image_size(2)
        current_points = nonzeros(points(i, j, :));
        num_points = length(current_points);
        if num_points > 0
            all_points(point_counter:point_counter + num_points - 1) = current_points;
            pixel_indices(point_counter:point_counter + num_points - 1, :) = repmat([i, j], num_points, 1);
            point_counter = point_counter + num_points;
        end
    end
end

% Perform k-means clustering on the flattened points
waitbar(3/4, f, "Clustering");
% left = 1;
% right = length(all_points);
first_el = squeeze(points(1, 1, :));
mins = first_el; % Contains the minimums for each group
groups = arrayfun(@(x) {x}, first_el(:));
clusters = zeros(length(all_points), 1);
clusters(1:length(first_el)) = 1:length(first_el);
curr_pixel = true(max_points, 1);
curr_dists = inf * ones(max_points, 1);
curr_pixel_x = 1;
curr_pixel_y = 1;
for k=(length(mins)+1):length(all_points) % for each peak after the first pixel
    i = pixel_indices(k, 1);
    j = pixel_indices(k, 2);
    if curr_pixel_x ~= i || curr_pixel_y ~= j
        curr_pixel_x = i;
        curr_pixel_y = j;
        curr_pixel = true(max_points, 1);
        curr_dists = inf * ones(max_points, 1);
    end
    dist = abs(mins - all_points(k));
    dist(dist >= curr_dists) = inf;

    % dist(~curr_pixel) = inf;
    [val, group] = min(dist);

    curr_pixel(group) = false;
    curr_dists(group) = val;
    mins(group) = all_points(k);
    groups{group} = [groups{group}, all_points(k)];
    clusters(k) = group;
end

last_el = squeeze(points(180, 180, :));
mins2 = last_el; % Contains the minimums for each group
groups2 = arrayfun(@(x) {x}, last_el(:));
clusters2 = zeros(length(all_points), 1);
clusters2((length(all_points)-length(last_el) + 1):length(all_points)) = 1:length(last_el);
curr_pixel = true(max_points, 1);
curr_pixel_x = size(points, 1);
curr_pixel_y = size(points, 2);
for k = (length(all_points) - length(mins2)):-1:1 % for each peak after the first pixel
    i = pixel_indices(k, 1);
    j = pixel_indices(k, 2);
    
    if curr_pixel_x ~= i || curr_pixel_y ~= j
        curr_pixel_x = i;
        curr_pixel_y = j;
        curr_pixel = true(max_points, 1);
    end
    
    % Find the closest element in mins2 to the current peak in all_points
    dist = abs(mins2 - all_points(k));
    dist(~curr_pixel) = inf; % Set distances of already used groups to infinity
    [~, group] = min(dist);
    
    % Update the group and mark it as used
    curr_pixel(group) = false;
    mins2(group) = all_points(k);
    groups2{group} = [groups2{group}, all_points(k)];
    clusters2(k) = group;
    
    if i == 90
        disp([abs(mins2 - all_points(k))', all_points(k)]);
        disp([curr_pixel_y, group, k, curr_pixel']);
    end
end

figure;
hold on;
colors = ['r', 'g', 'b', 'c', 'm', 'k'];
for k=1:length(all_points) % Horizontal line cut
    if pixel_indices(k, 1) == 90
        % disp([pixel_indices(k, 1), pixel_indices(k, 2), clusters(k)])
        scatter(all_points(k), pixel_indices(k, 2), colors(clusters(k)))
    end
end
set(gca, 'YDir','reverse')

figure;
hold on;
colors = ['r', 'g', 'b', 'c', 'm', 'k'];
for k=1:length(all_points) % Horizontal line cut
    if pixel_indices(k, 1) == 90
        % disp([pixel_indices(k, 2), clusters2(k)])
        % disp([pixel_indices(k, 1), pixel_indices(k, 2), clusters(k)])
        scatter(all_points(k), pixel_indices(k, 2), colors(clusters2(k)))
    end
end
set(gca, 'YDir','reverse')
disp([mins, mins2])

% seeds = zeros(max_points, 1);
% for k=1:max_points
%     seeds(k) = mean2(locs_array(:, :, k));
% end
% [cluster_indices, cluster_centers] = kmeans(all_points, max_points,"Distance","cityblock", "Start",seeds);

waitbar(7/8, f, "Assigning");
bucketed_points = zeros(size(points));
num_dupes = ones(max_points);
icache = 1;
jcache = 1;
% Assign clustered points back to their respective pixels and buckets
for k = 1:length(all_points)
    i = pixel_indices(k, 1);
    j = pixel_indices(k, 2);
    if i ~= icache || j ~= jcache % If change in pixel reset dupes
        icache = i;
        jcache = j;
        num_dupes = ones(max_points);
    end

    bucket_idx = cluster_indices(k);
    val = all_points(k);
    if bucketed_points(i, j, bucket_idx) ~= 0 % If cluster dupe
        % Store the average
        temp_total = val * num_dupes(bucket_idx) + bucketed_points(i, j, bucket_idx);
        num_dupes(bucket_idx) = num_dupes(bucket_idx) + 1;
        val = temp_total / num_dupes(bucket_idx);
    end
    % fails = false(size(cluster_centers));
    % while bucketed_points(i, j, bucket_idx) ~= 0
    %     fails(bucket_idx) = true;
    %     if all(fails)
    %         break
    %     end
    %     distances = abs(cluster_centers - all_points(k));
    %     distances(fails) = inf;
    %     [~, bucket_idx] = min(distances);
    % end
    bucketed_points(i, j, bucket_idx) = val;
end

mins = zeros(max_points, 1);
for k = 1:max_points
    data = squeeze(bucketed_points(:, :, k));
    data(data == 0) = NaN;
    mins(k) = min(data(:));    
end
[sorted_mins, indices] = sort(mins);

sorted_bucketed_points = zeros(image_size(1), image_size(2), max_points);
for k = 1:max_points
    sorted_bucketed_points(:, :, k) = bucketed_points(:, :, indices(k));
end

waitbar(7/8, f, "Plotting");
disp("old")
plot_point_pages(max_points, points)

disp("new")
plot_point_pages(max_points, sorted_bucketed_points)
close(f)
