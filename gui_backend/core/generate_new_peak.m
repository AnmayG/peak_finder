function [new_location, new_value, new_width, new_prom, error_coeff] = generate_new_peak(old_location, ...
    region_center, signal, raw, freq, old_locs, params_struct)
    % Generate a new peak location by getting all peaks in the area 
    % and choosing the one closest to the ideal
    % (after smoothing first, then without)
    tolerance = params_struct.tolerance;
    method = params_struct.generating_method;
    min_dist = params_struct.diff_peak_distance;
    error_coeff = 1;
    % region_center = 2 * shift - old_location;
    region_left = max([region_center - tolerance, freq(1)]);
    region_right = min([freq(end), region_center + tolerance]);
    region_indices = (freq >= region_left) & (freq <= region_right);
    signal_cropped = signal(region_indices);
    x_cropped = freq(region_indices);

    if isempty(x_cropped) || length(x_cropped) <= 3
        new_location = -1;
        new_value = -1;
        new_width = -1;
        new_prom = -1;
        return
    end
    
    if method == 7
        new_location = -1;
        new_value = -1;
        new_width = -1;
        new_prom = -1;
        return
    end

    threshold = -min(-signal_cropped);
    baseline = mean(signal(1:10));

    if method == 1 % Tallest, if split average
        [vals, locs, widths, proms] = ...
            findpeaks_custom(-signal_cropped + threshold, x_cropped, ...
                            'SortStr','descend', ...
                            'MinPeakDistance', min_dist, ...
                            'WidthReference','halfheight', ...
                            'NPeaks',2);
        vals = -vals + threshold;
        % [~, prom_index] = maxk(vals, 1);
        % vals = vals(prom_index);
        % locs = locs(prom_index);
        % widths = widths(prom_index);
        % proms = proms(prom_index);

        % if numel(locs) < num_max_peaks
        % Search for positive contrast and include that as well
            upper_peak_threshold = prctile(signal_cropped, 90); % Top X%
            contrast = (max(signal_cropped) - baseline) / baseline;
            if contrast > 0.002 % Positive contrast detected
                [vals2, locs_peaks2, widths2, prominences2] = findpeaks_custom(signal_cropped - upper_peak_threshold, x_cropped, ...
                                    'SortStr','descend',...
                                    'WidthReference','halfheight', ...
                                    'NPeaks', 1);
                vals2 = vals2 + upper_peak_threshold;
                
                vals = [vals; vals2];
                locs = [locs, locs_peaks2];
                widths = [widths, widths2];
                proms = [proms; prominences2];
            end
        % end
    
        % Duplicates from already found peaks are removed
        [~, duplicate_index] = intersect(locs, old_locs);
        locs(duplicate_index) = [];
        vals(duplicate_index) = [];
        widths(duplicate_index) = [];
        proms(duplicate_index) = [];
        failure = isempty(locs);
    
        % Actual peak is split or flattened
        if numel(locs) == 2
            contrast = abs(vals - baseline) ./ baseline;
            if abs(vals(1) - vals(2)) < 2
                if contrast(1) > 0.3e-3 && contrast(2) > 0.3e-3
                    failure = true;
                end
            else
                [~, prom_index] = mink(vals, 1);
                vals = vals(prom_index);
                locs = locs(prom_index);
                widths = widths(prom_index);
                proms = proms(prom_index);
            end
        end

        % Consider average of two peaks
        if failure
            error_coeff = 2;
            y_data = -signal_cropped + threshold;
            [vals, locs, widths, proms] = ...
                findpeaks_custom(y_data, x_cropped, ...
                                'SortStr','descend',...
                                'WidthReference','halfheight', ...
                                'NPeaks',2);
            vals = -vals + threshold; % flip it back
            
            duplicate_index = find(locs == old_location);
            locs(duplicate_index) = [];
            vals(duplicate_index) = [];
            widths(duplicate_index) = [];
            proms(duplicate_index) = [];
    
            locs = mean(locs);
            vals = mean(vals);
            widths = mean(widths);
            proms = mean(proms);
        end
    elseif method == 2 % Highest prominence
        [vals, locs, widths, proms] = ...
            findpeaks_custom(-signal_cropped + threshold, x_cropped, ...
                            'SortStr','descend',...
                            'MinPeakDistance', min_dist, ...
                            'WidthReference','halfheight');
        vals = -vals + threshold; % flip it back
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
    elseif method == 3 % Highest goodness
        [vals, locs, widths, proms] = ...
            findpeaks_custom(-signal_cropped + threshold, x_cropped, ...
                            'SortStr','descend',...
                            'MinPeakDistance', min_dist, ...
                            'WidthReference','halfheight');
        vals = -vals + threshold; % flip it back
        [~, duplicate_index] = intersect(locs, old_locs);
        locs(duplicate_index) = [];
        vals(duplicate_index) = [];
        widths(duplicate_index) = [];
        proms(duplicate_index) = [];
    elseif method == 4 % Immediately average
        y_data = -signal_cropped + threshold;
        [vals, locs, widths, proms] = ...
                findpeaks_custom(y_data, x_cropped, ...
                                'SortStr','descend',...
                                'MinPeakDistance', min_dist, ...
                                'WidthReference','halfheight', ...
                                'NPeaks',2);
        vals = -vals + threshold; % flip it back
        
        duplicate_index = find(locs == old_location);
        locs(duplicate_index) = [];
        vals(duplicate_index) = [];
        widths(duplicate_index) = [];
        proms(duplicate_index) = [];

        locs = mean(locs);
        vals = mean(vals);
        widths = mean(widths);
        proms = mean(proms);
    elseif method == 5 % Immediately Average of two proms
        y_data = -signal_cropped + threshold;
        [vals, locs, widths, proms] = ...
                findpeaks_custom(y_data, x_cropped, ...
                                'SortStr','descend',...
                                'MinPeakDistance', min_dist, ...
                                'WidthReference','halfheight');
        vals = -vals + threshold; % flip it back
        [~, prom_index] = maxk(proms, 2);
        vals = vals(prom_index);
        locs = locs(prom_index);
        widths = widths(prom_index);
        proms = proms(prom_index);
        
        duplicate_index = find(locs == old_location);
        locs(duplicate_index) = [];
        vals(duplicate_index) = [];
        widths(duplicate_index) = [];
        proms(duplicate_index) = [];

        locs = mean(locs);
        vals = mean(vals);
        widths = mean(widths);
        proms = mean(proms);
    elseif method == 6 % Tallest
        [vals, locs, widths, proms] = ...
            findpeaks_custom(-signal_cropped + threshold, x_cropped, ...
                            'SortStr','descend',...
                            'MinPeakDistance', min_dist, ...
                            'WidthReference','halfheight', ...
                            'NPeaks', 1);
        vals = -vals + threshold;
        % [~, prom_index] = maxk(vals, 1);
        % vals = vals(prom_index);
        % locs = locs(prom_index);
        % widths = widths(prom_index);
        % proms = proms(prom_index);
    
        % Duplicates from already found peaks are removed
        [~, duplicate_index] = intersect(locs, old_locs);
        locs(duplicate_index) = [];
        vals(duplicate_index) = [];
        widths(duplicate_index) = [];
        proms(duplicate_index) = [];
    end

    % Zoning constraint 1
    tmp_peaks_info = struct('locs', locs, 'vals', vals, 'widths', widths, 'proms', proms);
    rm_idx = filter_peaks(tmp_peaks_info, params_struct, old_location);
    locs(rm_idx) = [];

    % Ultimate failure case: consider the raw
    % Consider lightly smoothed as per Fu's code
    if isempty(locs)
        error_coeff = 3;
        raw_cropped = raw(region_indices);
        z = smooth(raw_cropped, 'sgolay', 2);
        [vals, locs, widths, proms] = ...
            findpeaks_custom(-z + threshold, x_cropped, ...
                            'SortStr','descend',...
                            'MinPeakDistance', min_dist, ...
                            'WidthReference','halfheight');
        vals = -vals + threshold;
        
        duplicate_index = find(locs == old_location);
        locs(duplicate_index) = [];
        vals(duplicate_index) = [];
        widths(duplicate_index) = [];
        proms(duplicate_index) = [];

        % Consider raw
        if isempty(locs)
            error_coeff = 4;
            [vals, locs, widths, proms] = ...
                findpeaks_custom(-raw_cropped + threshold, x_cropped, ...
                                'SortStr','descend',...
                                'WidthReference','halfheight');
            vals = -vals + threshold;
            duplicate_index = find(locs == old_location);
            locs(duplicate_index) = [];
            vals(duplicate_index) = [];
            widths(duplicate_index) = [];
            proms(duplicate_index) = [];
        end
    end
    
    % Zoning constraint 2
    tmp_peaks_info2 = struct("locs", locs, "vals", vals, "widths", widths, "proms", proms);
    rm_idx2 = filter_peaks(tmp_peaks_info2, params_struct, old_location);
    locs(rm_idx2) = [];
    % If still nothing, return none
    if isempty(locs)
        error_coeff = 4;
        new_location = -1;
        new_value = -1;
        new_width = -1;
        new_prom = -1;
    else
        % If a list exists, return the one with the highest goodness
        peaks_goodness = -(locs - region_center).^2;
        [~, new_ind] = max(peaks_goodness);
        new_location = locs(new_ind);
        new_value = vals(new_ind);
        new_width = widths(new_ind);
        new_prom = proms(new_ind);
    end
end
