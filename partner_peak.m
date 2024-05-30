function [locs, vals, widths, index_pairs, error_coeff]...
          = partner_peak(locs, ...,
                         split, ...
                         tolerance, ...
                         signal, ...
                         threshold, ...
                         freq, ...
                         vals, widths, proms, display, raw)
    if nargin < 10
        display = false;
    end
    % Check to see if a real partner peak exists by looking in the possible
    % physical interval
    [vals, sorted_indices] = sort(vals); % Sort peaks from left to right
    sorted_indices = sorted_indices';
    locs = locs(sorted_indices); % Rearrange vals and widths accordingly
    widths = widths(sorted_indices);
    proms = proms(sorted_indices);
    goodness_matrix = zeros(length(sorted_indices), length(sorted_indices)); % Adjacency matrix between peaks for matching
    for i = sorted_indices % conventional indexes still work!
        disp(i)
        peak = locs(i);
        peak_goodness = partner_peak_vector(peak, locs, split, tolerance);
        goodness_matrix(i, :) = peak_goodness;
    end

    % we love chatgpt
    [num_rows, num_cols] = size(goodness_matrix);
    index_pairs = zeros(num_rows, 2);
    
    values_list = zeros(num_rows * num_cols, 3);
    ind = 1;
    for i = 1:num_rows
        for j = 1:num_cols
            new_row = [i, j, goodness_matrix(i, j)];
            values_list(ind, :) = new_row;
            ind = ind + 1;
        end
    end
    values_list = sortrows(values_list, -3);
    
    assigned_rows = false(1, num_rows);
    assigned_cols = false(1, num_cols);
    assigned_count = 0;
    for k = 1:size(values_list, 1)
        row = values_list(k, 1);
        col = values_list(k, 2);
        
        if ~assigned_rows(row) && ~assigned_cols(col)
            index_pairs(row, :) = [row, col];
            assigned_rows(row) = true;
            assigned_cols(col) = true;
            assigned_count = assigned_count + 1;
        end
        
        % Stop if all rows have been assigned
        if assigned_count == num_rows
            break;
        end
    end

    error_coeff = 4;
    [locs, vals, widths, index_pairs, error_coeff] = ...
        add_new_peak(locs', vals, widths', proms', index_pairs, ...
        split, tolerance, signal, raw, freq);

    if display
        figure
        hold on
        yline(threshold, 'Color', 'r')
        plot(freq, signal, 'Color', 'b')
        plot(locs, vals, 'rs', 'MarkerFaceColor', 'k')
        colors = ['r', 'g', 'b', 'c', 'm', 'y'];

        for i = 1:size(index_pairs, 1) % For each pair
            if i <= ceil(size(index_pairs, 1) / 2)
                peak = locs(index_pairs(i, 1)); % Get the first peak and second peak
                peak_partner = locs(index_pairs(i, 2));
                init_height = vals(index_pairs(i, 1));
                second_height = vals(index_pairs(i, 2));
                color = colors(mod(i - 1, 6) + 1); % Give them a color
                plot(peak, init_height, 'rv', 'MarkerFaceColor', color)
                plot(peak_partner, second_height, 'rv', 'MarkerFaceColor', color)
                text(peak_partner + 0.02e9,second_height,num2str(i))
                text(peak + 0.02e9,init_height,num2str(i))
    
                xline(split, 'LineWidth', 2)
                if i == 1
                    center = 2 * split - peak;
                    xline(peak)
                    xline(center)
                    xregion(center - tolerance, center + tolerance, 'FaceColor', color)
                    % center = 2 * split - peak_partner;
                    % xline(center)
                    % xregion(center - tolerance, center + tolerance, 'FaceColor', color)
                end
            end
        end

        % title(num2str(index_pairs(:, 2))')
        hold off
    end
end

% peak_loc is the initial peak that is looking for a partner
% locs is list of all peak locations
% split is the dividing point that a peak should be centered at
% tolerance is the tolerance from the ideal (physical) point
% Must guarantee that peak_loc < split otherwise it'll pick itself
% Returns an array of "goodness", a value from -Inf - 0 representing how
% "good" a peak is compared to its physical ideal
% Higher "goodness" means that it's a better ideal
% This allows us to determine if an imaginary peak is necessary by the
% pigeonhole principle; if only the best pairs exist, then any peaks left 
% unpaired likely have an undiscovered good partner (or they're noise
% peaks and thus should be ignored, but these would ideally come with a warning of some kind)
function [peaks_goodness] = partner_peak_vector(peak_loc, ...
                                                locs, ...
                                                split, ...
                                                tolerance)
    % Calculate the range of values that a peak can physically exist
    phantom_peak_center = 2 * split - peak_loc; % ideal partner location
    phantom_peak_left = phantom_peak_center - tolerance;
    phantom_peak_right = phantom_peak_center + tolerance;
    % Iterate through all of the peaks and find "goodness" for each one
    peaks_goodness = -Inf * ones(size(locs));
    for i=1:numel(locs)
        peak = locs(i);
        if peak ~= peak_loc % Given peak is on left of split, partner is on right
            if peak > phantom_peak_left && peak < phantom_peak_right
                % Inside possible interval, calculate goodness of each peak
                peaks_goodness(i) = -(peak - phantom_peak_center)^2;
            end
        end
    end
end

function [updated_locations, updated_values, updated_widths, updated_pairs, error_coeff] = ...
    add_new_peak(locations, values, widths, proms, pairs, split, tolerance, signal, raw, x)
    % Initialize updated arrays with the original lists
    updated_locations = locations;
    updated_values = values;
    updated_widths = widths;
    updated_pairs = pairs; % Copy of original pairs
    error_coeff = 0;
        
    % Identify self-pairs and generate new locations, values, and widths
    mark_for_deletion = [];
    for i = 1:size(pairs, 1)
        if pairs(i, 1) == pairs(i, 2)
            % Self-pair found, generate new location, value, and width
            [new_location, new_value, new_width, error_coeff] = ...
                generate_new_peak(locations(pairs(i, 1)), split, ...
                tolerance, signal, raw, x, updated_locations, proms(pairs(i, 1)));

            if new_location == -1
                % No valid pair found, peak should be deleted
                mark_for_deletion = [mark_for_deletion, locations(pairs(i, 1))];
                continue;
            end

            % Determine the new index for the new location
            new_index = find_insertion_index(updated_locations, new_location);
            
            % Update all pairs ahead of the insert
            for j=1:size(updated_pairs, 1)
                if updated_pairs(j, 1) >= new_index
                    updated_pairs(j, 1) = updated_pairs(j, 1) + 1;
                end
                if updated_pairs(j, 2) >= new_index
                    updated_pairs(j, 2) = updated_pairs(j, 2) + 1;
                end
            end

            % Remember the index of the self-pair to replace
            replace_index = find(updated_locations == locations(pairs(i, 1)));
            replace_index = replace_index(1);
            updated_pairs(replace_index, 2) = new_index;

            updated_locations = insert_sorted(updated_locations, new_location, new_index);
            updated_values = insert_sorted(updated_values, new_value, new_index);
            updated_widths = insert_sorted(updated_widths, new_width, new_index);
            if replace_index >= new_index % Check to see if inserting messed anything up
                replace_index = replace_index + 1;
            end

            % Add in the backwards edge as well
            new_row = [new_index, replace_index];
            updated_pairs = [updated_pairs(1:new_index - 1, :); new_row; updated_pairs(new_index:end, :)];
        end
    end

    for bad_peak=mark_for_deletion
        i = find(updated_locations == bad_peak);
        updated_locations(i) = [];
        updated_values(i) = [];
        updated_widths(i) = [];
        updated_pairs(i, :) = [];
        % Modify after delete,
        for j=1:size(updated_pairs, 1)
                if updated_pairs(j, 1) >= i
                    updated_pairs(j, 1) = updated_pairs(j, 1) - 1;
                end
                if updated_pairs(j, 2) >= i
                    updated_pairs(j, 2) = updated_pairs(j, 2) - 1;
                end
        end
    end
end

function [new_location, new_value, new_width, error_coeff] = generate_new_peak(old_location, ...
    split, tolerance, signal, raw, x, old_locs, old_prom)
    % Generate a new peak location by getting all peaks in the area 
    % and choosing the one closest to the ideal
    % (after smoothing first, then without)
    % new_location = old_location + 0.01e9; % Adjust the criteria as needed
    % new_value = 2.7e4;
    % new_width = 1e9;

    error_coeff = 1;
    region_center = 2 * split - old_location;
    region_left = max([region_center - tolerance, x(1)]);
    region_right = min([x(end), region_center + tolerance]);
    region_indices = (x >= region_left) & (x <= region_right);
    signal_cropped = signal(region_indices);
    x_cropped = x(region_indices);
    if isempty(x_cropped) || length(x_cropped) <= 3
        new_location = -1;
        new_value = -1;
        new_width = -1;
        return
    end
    threshold = -min(-signal_cropped);
    
    [vals, locs, widths, proms] = ...
        findpeaks(-signal_cropped + threshold, x_cropped, ...
                        'SortStr','descend',...
                        'WidthReference','halfheight');
    vals = -vals + threshold;
    [~, prom_index] = maxk(proms, 1);
    vals = vals(prom_index);
    locs = locs(prom_index);
    widths = widths(prom_index);
    proms = proms(prom_index);
    [~, duplicate_index] = intersect(locs, old_locs);
    locs(duplicate_index) = [];
    vals(duplicate_index) = [];
    widths(duplicate_index) = [];
    proms(duplicate_index) = [];
    % Consider lightly smoothed as per Fu's code
    if isempty(locs)
        error_coeff = 2;
        raw_cropped = raw(region_indices);
        z = smooth(raw_cropped, 'sgolay', 2);
        [vals, locs, widths, proms] = ...
            findpeaks(-z + threshold, x_cropped, ...
                            'SortStr','descend',...
                            'WidthReference','halfheight');
        vals = -vals + threshold;
        duplicate_index = find(locs == old_location);
        locs(duplicate_index) = [];
        vals(duplicate_index) = [];
        widths(duplicate_index) = [];
        proms(duplicate_index) = [];

        % Consider raw
        if isempty(locs)
            error_coeff = 3;
            [vals, locs, widths, proms] = ...
                findpeaks(-raw_cropped + threshold, x_cropped, ...
                                'SortStr','descend',...
                                'WidthReference','halfheight');
            vals = -vals + threshold;
            duplicate_index = find(locs == old_location);
            locs(duplicate_index) = [];
            vals(duplicate_index) = [];
            widths(duplicate_index) = [];
            proms(duplicate_index) = [];

            % Just keep on widening the threshold
            % while isempty(locs)
            %     threshold = threshold * 5;
            %     [vals, locs, widths] = ...
            %         findpeaks(-raw_cropped + threshold, x_cropped, ...
            %                         'SortStr','descend',...
            %                         'WidthReference','halfheight');
            %     vals = -vals + threshold;
            %     duplicate_index = find(locs == old_location);
            %     locs(duplicate_index) = [];
            %     vals(duplicate_index) = [];
            %     widths(duplicate_index) = [];
            % end
        end
        % figure
        % hold on
        % plot(x, raw, x_cropped, raw_cropped)
        % plot(locs, vals, 'rv')
        % hold off
    end
    if isempty(locs)
        error_coeff = 4;
        new_location = -1;
        new_value = -1;
        new_width = -1;
    else
        peaks_goodness = -(locs - region_center).^2;
        [~, new_ind] = max(peaks_goodness);
        new_location = locs(new_ind);
        new_value = vals(new_ind);
        new_width = widths(new_ind);
    end

    % figure
    % hold on
    % plot(x, raw)
    % plot(locs(new_ind), vals(new_ind), 'rv')
    % hold off
end

function sorted_array = insert_sorted(array, new_element, idx)
    % Insert the new element while keeping the array sorted
    if isempty(idx)
        sorted_array = [array; new_element]; % Add to the end
    else
        sorted_array = [array(1:idx-1); new_element; array(idx:end)];
    end
end

function idx = find_insertion_index(locations, new_location)
    % Find the correct index to insert the new location to keep the list sorted
    idx = find(locations >= new_location, 1);
    if isempty(idx)
        idx = size(locations, 1) + 1; % Add to the end
    end
end

% load('2p03GPa_14p9k_23db_39G.mat');
% 
% x_grid_width = 64;
% y_grid_width = 64;
% for x = randi([0, x_grid_width - 1], 1, 256 / x_grid_width) + (1:x_grid_width:256)
%     for y = randi([0, y_grid_width - 1], 1, 256 / y_grid_width) + (1:y_grid_width:256)
%         [vals, locs, widths, proms, threshold, signal] = find_peaks_at_point(x, y, gWide);
%         split = 2.91e9;
%         tolerance = 0.05e9;
%         [index_pairs, sorted_locs, sorted_indices] = partner_peak2(locs, split, tolerance, false, signal, threshold, gWide, vals);
%     end
% end
