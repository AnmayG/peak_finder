function [vals, locs, widths, proms, peak_threshold, error] = ...
    find_peaks_at_point(z, freq, display, num_max_peaks, peak_perc_threshold, ...
    diff_peak_distance)
        if nargin < 3 % Don't need to include display if not necessary
            display = false;
            num_max_peaks = 2;
            % Previously 32%
            peak_perc_threshold = 5;
            diff_peak_distance = 0.02e9;
        end
        error = 0;
        % Chop off the ends to minimize chance of badness happening - plus
        % peaks shouldn't be there anyways
        x_bounds = (freq > freq(1) + diff_peak_distance) & (freq < freq(end) - diff_peak_distance);
        
        peak_threshold = prctile(z, peak_perc_threshold);
        [vals, locs, widths, proms] = findpeaks_custom(max(-z(x_bounds) + peak_threshold, 0), freq(x_bounds), ...
                            'SortStr','descend',...
                            'WidthReference','halfheight', ...
                            'NPeaks', num_max_peaks, ...,
                            'MinPeakDistance', diff_peak_distance);
        vals = -vals + peak_threshold;
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

        if display
            figure
            hold on
            yline(peak_threshold, 'Color', 'r')
            plot(freq, z, 'Color', 'b')
            plot(locs, vals, 'rv', 'MarkerFaceColor', 'g');
            title("Point at " + x + ", " + y)
            disp(x + " " + y + " completed");
        end
end
