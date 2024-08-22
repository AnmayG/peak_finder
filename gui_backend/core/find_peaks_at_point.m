function peaks_info = ...
    find_peaks_at_point(signal, freq, display, params_struct)
        if nargin < 3 % Don't need to include display if not necessary
            display = false;
            params_struct = struct();
            params_struct.num_max_peaks = 2;
            % Previously 32%
            params_struct.peak_perc_threshold = 5;
            params_struct.diff_peak_distance = 0.02e9;
        end
        error = 0;
        % Chop off the ends to minimize chance of badness happening - plus
        % peaks shouldn't be there anyways
        % x_bounds = (freq > freq(1) + params_struct.diff_peak_distance) & (freq < freq(end) - params_struct.diff_peak_distance);
        x_bounds = (freq > freq(1)) & (freq < freq(end));
        
        if params_struct.normalize
            peak_threshold = mean(signal) - params_struct.peak_perc_threshold * std(signal);
        else
            peak_threshold = prctile(signal, params_struct.peak_perc_threshold);
        end
        [vals, locs, widths, proms] = findpeaks_custom(max(-signal(x_bounds) + peak_threshold, 0), freq(x_bounds), ...
                            'SortStr','descend',...
                            'WidthReference','halfheight', ...
                            'NPeaks', params_struct.num_max_peaks, ...,
                            'MinPeakDistance', params_struct.diff_peak_distance);
        vals = -vals + peak_threshold;
        % Somehow these get transposed idk why
        vals = vals';
        proms = proms';
        % Remove all peaks with abnormally low contrast
        % temp_baseline = mean(z(1:10)); % average of first 10 is a bad baseline, replace with something smarter later
        % contrasts = (vals - temp_baseline) ./ temp_baseline;
        % rm_locs = contrasts > -0.002;
        % if sum(rm_locs) < numel(locs)
        %     error = 1;
        %     vals(rm_locs) = [];
        %     locs(rm_locs) = [];
        %     widths(rm_locs) = [];
        %     proms(rm_locs) = [];
        % end

        % Sort peaks from left to right
        [~, sorted_indices] = sort(locs);
        sorted_indices = sorted_indices';
        locs = locs(sorted_indices)';
        vals = vals(sorted_indices)'; % Rearrange accordingly
        widths = widths(sorted_indices)';
        proms = proms(sorted_indices)';
        
        peaks_info = struct();
        peaks_info.vals = vals;
        peaks_info.locs = locs;
        peaks_info.widths = widths;
        peaks_info.proms = proms;
        peaks_info.threshold = peak_threshold;
        peaks_info.error = error;

        if display
            figure
            hold on
            yline(peak_threshold, 'Color', 'r')
            plot(freq, signal, 'Color', 'b')
            plot(locs, vals, 'rv', 'MarkerFaceColor', 'g');
            title("Point at " + x + ", " + y)
            disp(x + " " + y + " completed");
        end
end
