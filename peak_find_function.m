function seed = peak_find_function(split, tolerance, x, y, data, freq, num_max_peaks)
        if nargin < 7
            num_max_peaks = 2;
            peak_perc_threshold = 5;
            diff_peak_distance = 0.002;
            smooth_span = 41;
            smooth_degree = 7;
        end
        [vals, locs, widths, proms, threshold, signal, raw] = ...
            find_peaks_at_point(y, x, data, freq, false, num_max_peaks, ...
            peak_perc_threshold, diff_peak_distance, smooth_span, smooth_degree);
        [sorted_locs, sorted_vals, sorted_widths, pairs, error_coeff] = ...
            partner_peak(locs, split, tolerance, signal, threshold, freq, vals, widths, proms, 0, false, raw);
        baseline = find_baseline(y, x, data, freq);
        seed = peaks_to_seed(sorted_locs, sorted_vals, sorted_widths, pairs, baseline, error_coeff);
end
