function [new_location, new_value, new_width, error_coeff] = generate_new_peak(old_location, ...
    split, tolerance, signal, raw, freq, old_locs, method)
    % Generate a new peak location by getting all peaks in the area 
    % and choosing the one closest to the ideal
    % (after smoothing first, then without)
    % new_location = old_location + 0.01e9; % Adjust the criteria as needed
    % new_value = 2.7e4;
    % new_width = 1e9;

    error_coeff = 1;
    region_center = 2 * split - old_location;
    region_left = max([region_center - tolerance, freq(1)]);
    region_right = min([freq(end), region_center + tolerance]);
    region_indices = (freq >= region_left) & (freq <= region_right);
    signal_cropped = signal(region_indices);
    x_cropped = freq(region_indices);

    if isempty(x_cropped) || length(x_cropped) <= 3
        new_location = -1;
        new_value = -1;
        new_width = -1;
        return
    end

    threshold = -min(-signal_cropped);
    baseline = mean(signal(1:10));

    if method == 1 % Tallest, if split average
        [vals, locs, widths, proms] = ...
            findpeaks(-signal_cropped + threshold, x_cropped, ...
                            'SortStr','descend',...
                            'WidthReference','halfheight', ...
                            'NPeaks',2);
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
        failure = isempty(locs);
    
        % Actual peak is split or flattened
        if numel(locs) == 2
            contrast = abs(vals - baseline) ./ baseline;
            if abs(vals(1) - vals(2)) < 2
                if contrast(1) > 0.003 && contrast(2) > 0.5e-4
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
                findpeaks(y_data, x_cropped, ...
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
            findpeaks(-signal_cropped + threshold, x_cropped, ...
                            'SortStr','descend',...
                            'WidthReference','halfheight');
        vals = -vals + threshold; % flip it back
        [~, prom_index] = maxk(proms, 1);
        vals = vals(prom_index);
        locs = locs(prom_index);
        widths = widths(prom_index);
        proms = proms(prom_index);
    elseif method == 3 % Highest goodness
        [vals, locs, widths, ~] = ...
            findpeaks(-signal_cropped + threshold, x_cropped, ...
                            'SortStr','descend',...
                            'WidthReference','halfheight');
        vals = -vals + threshold; % flip it back
    elseif method == 4 % Immediately average
        y_data = -signal_cropped + threshold;
        [vals, locs, widths, proms] = ...
                findpeaks(y_data, x_cropped, ...
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
    elseif method == 5 % Immediately Average of two proms
        y_data = -signal_cropped + threshold;
        [vals, locs, widths, proms] = ...
                findpeaks(y_data, x_cropped, ...
                                'SortStr','descend',...
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
            findpeaks(-signal_cropped + threshold, x_cropped, ...
                            'SortStr','descend',...
                            'WidthReference','halfheight', ...
                            'NPeaks',1);
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
        failure = isempty(locs);
    end

    % Ultimate failure case: consider the raw
    % Consider lightly smoothed as per Fu's code
    if isempty(locs)
        error_coeff = 3;
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
            error_coeff = 4;
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
        end
    end

    % If still nothing, return none
    if isempty(locs)
        error_coeff = 4;
        new_location = -1;
        new_value = -1;
        new_width = -1;
    else
        % If a list exists, return the one with the highest goodness
        peaks_goodness = -(locs - region_center).^2;
        [~, new_ind] = max(peaks_goodness);
        new_location = locs(new_ind);
        new_value = vals(new_ind);
        new_width = widths(new_ind);
    end
end
