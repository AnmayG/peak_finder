function [vals, locs_peaks, widths, prominences, peak_threshold, z, raw] = ...
    find_peaks_at_point(y, x, data, freq, display, num_max_peaks, peak_perc_threshold, ...
    diff_peak_distance, smooth_span, smooth_degree)
        if nargin < 5 % Don't need to include display if not necessary
            display = false;
            num_max_peaks = 2;
            % Previously 32%
            peak_perc_threshold = 5;
            diff_peak_distance = 0.02e9;
            smooth_span = 41;
            smooth_degree = 7;
        end

        raw = squeeze(data(x, y, :));
        % z = sgolayfilt(z, 7, 41);
        z = smooth(raw, smooth_span, 'sgolay', smooth_degree);
        peak_threshold = prctile(z, peak_perc_threshold);
        % offset = min(-z); % Lowest point in flipped version
        [vals, locs_peaks, widths, prominences] = findpeaks(-z+peak_threshold, freq, ...
                            'SortStr','descend',...
                            'WidthReference','halfheight', ...
                            'NPeaks', num_max_peaks, ...,
                            'MinPeakDistance', diff_peak_distance);
        vals = -vals + peak_threshold;
        if display
            figure
            hold on
            yline(peak_threshold, 'Color', 'r')
            plot(gWide.SweepParam, z, 'Color', 'b')
            plot(locs_peaks, vals, 'rv', 'MarkerFaceColor', 'g');
            title("Point at " + x + ", " + y)
            disp(x + " " + y + " completed");
        end
end
