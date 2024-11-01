function new_peaks_info ...
          = partner_peak(signal, ...
                         raw, ...
                         freq, ...
                         peaks_info, ...
                         params_struct, ...
                         display)
    if nargin < 7
        display = false;
    end

    % Extract info from structs
    centers = peaks_info.centers;
    locs = peaks_info.locs;
    tolerance = params_struct.tolerance;
    if(isfield(params_struct, "zones_x"))
        zones_x = params_struct.zones_x;
        zones_y = params_struct.zones_y;
        if(isempty(zones_x) || isempty(zones_y))
            zones_x = [0, Inf];
            zones_y = [0, Inf];
        end
        zoning_method = params_struct.zoning_method;
    else
        zones_x = [0, Inf];
        zones_y = [0, Inf];
        zoning_method = 0;
    end

    % Calculate goodness for partnering
    goodness_matrix = zeros(length(locs), length(locs)); % Adjacency matrix between peaks for matching
    for i = 1:length(locs)
        vec = partner_peak_vector(locs(i), locs, centers(i), tolerance);
        goodness_matrix(i, :) = vec;
    end

    % Zoning filters - If shifting/splitting based we don't want to pair
    % i.e. if we want shiftings from 3e7 to 5e7 and 2e9 to 3e9 we want to
    % keep the peaks but we don't want to partner them with one another
    for i = 1:length(locs)
        for j = 1:length(locs)
            peak1 = locs(i);
            peak2 = locs(j);
            if zoning_method == 4 || zoning_method == 7
                % Shifting
                peak_x = abs(peak1 + peak2) / 2;
                if ~is_valid_peak(peak_x, zones_x)
                    goodness_matrix(i, j) = -Inf;
                end
            end
            if zoning_method == 5 || zoning_method == 7
                % Splitting
                peak_x = abs(peak1 - peak2) / 2;
                if ~is_valid_peak(peak_x, zones_y)
                    goodness_matrix(i, j) = -Inf;
                end
            end
        end
    end

    % we love chatgpt
    [num_rows, num_cols] = size(goodness_matrix);
    index_pairs = zeros(num_rows, 2);
    
    % Transform the goodness matrix into a nice row-by-row format
    [row_indices, col_indices] = ndgrid(1:num_rows, 1:num_cols);
    values_list = [row_indices(:), col_indices(:), goodness_matrix(:)];

    % Sort them so that the best ones are higher
    values_list = sortrows(values_list, -3);
    
    % Assign rows by height
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
    peaks_info.pairs = index_pairs;
    new_peaks_info = add_new_peak(peaks_info, signal, raw, freq, params_struct);

    if display
        figure
        hold on
        yline(peaks_info.threshold, 'Color', 'r')
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
    
                if i == 1
                    center = centers(index_pairs(i, 1));
                    xline(peak)
                    xline(center)
                    xregion(center - tolerance, center + tolerance, 'FaceColor', color)
                end
            end
        end
        hold off
    end
end

function result = is_valid_peak(coord, valid_locs)
    result = false;
    num_intervals = length(valid_locs) - 1;
    for i = 1:2:num_intervals
        if coord >= valid_locs(i) && coord <= valid_locs(i+1)
            result = true;
            return;
        end
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
                                                phantom_peak_center, ...
                                                tolerance)
    % Calculate the range of values that a peak can physically exist
    % phantom_peak_center = 2 * reflectio - peak_loc; % ideal partner location
    phantom_peak_left = phantom_peak_center - tolerance;
    phantom_peak_right = phantom_peak_center + tolerance;
    % Iterate through all of the peaks and find "goodness" for each one
    peaks_goodness = -Inf * ones(size(locs));
    % for i=1:numel(locs)
    %     peak = locs(i);
    %     if peak ~= peak_loc
    %         if peak > phantom_peak_left && peak < phantom_peak_right
    %             % Inside possible interval, calculate goodness of each peak
    %             peaks_goodness(i) = -(peak - phantom_peak_center)^2;
    %         end
    %     end
    % end
    
    % Find the indices of peaks that are not the original and within the possible interval
    valid_indices = (locs ~= peak_loc) & (locs > phantom_peak_left) & (locs < phantom_peak_right);
    peaks_goodness(valid_indices) = -(locs(valid_indices) - phantom_peak_center).^2;
end

function new_peaks_info = add_new_peak(peaks_info, signal, raw, freq, params_struct)
    % Initialize updated arrays with the original lists
    locations = peaks_info.locs;
    updated_locations = peaks_info.locs;
    updated_values = peaks_info.vals;
    updated_widths = peaks_info.widths;
    updated_proms = peaks_info.proms;
    pairs = peaks_info.pairs;
    updated_pairs = peaks_info.pairs; % Copy of original pairs
    centers = peaks_info.centers;
    error_coeff = 0;
        
    % Identify self-pairs and generate new locations, values, and widths
    mark_for_deletion = [];
    for i = 1:size(pairs, 1)
        if pairs(i, 1) == pairs(i, 2)
            % Self-pair found, generate new location, value, and width

            % ideal_center = 2 * shift - locations(pairs(i, 1));
            ideal_center = centers(pairs(i, 1));

            [new_location, new_value, new_width, new_prom, error_coeff] = ...
                generate_new_peak(locations(pairs(i, 1)), ideal_center, ...
                signal, raw, freq, updated_locations, params_struct);

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
            updated_proms = insert_sorted(updated_proms, new_prom, new_index);
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
        updated_proms(i) = [];
        updated_pairs(i, :) = [];
        % Modify after delete
        for j=1:size(updated_pairs, 1)
                if updated_pairs(j, 1) >= i
                    updated_pairs(j, 1) = updated_pairs(j, 1) - 1;
                end
                if updated_pairs(j, 2) >= i
                    updated_pairs(j, 2) = updated_pairs(j, 2) - 1;
                end
        end
    end

    % Output packaging
    new_peaks_info = struct();
    new_peaks_info.locs = updated_locations;
    new_peaks_info.vals = updated_values;
    new_peaks_info.widths = updated_widths;
    new_peaks_info.proms = updated_proms;
    new_peaks_info.pairs = updated_pairs;
    new_peaks_info.error = error_coeff;
end

function sorted_array = insert_sorted(array, new_element, idx)
    % Insert the new element while keeping the array sorted    
    if isempty(idx)
        sorted_array = [array, new_element]; % Add to the end
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