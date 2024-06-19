function seed = peak_find_function(use_splitting, split, tolerance, x, y, data, freq, num_max_peaks, peak_perc_threshold, ...
    diff_peak_distance, smooth_span, smooth_degree, method, seed_method)
        if nargin < 8
            num_max_peaks = 2;
            peak_perc_threshold = 5;
            diff_peak_distance = 0.002;
            smooth_span = 41;
            smooth_degree = 7;
        end

        [vals, locs, widths, proms, threshold, signal, raw, error] = find_peaks_at_point(y, x, data, freq, ...
            false, num_max_peaks, peak_perc_threshold, diff_peak_distance, ...
            smooth_span, smooth_degree);
        
        centers = get_centers(locs, freq, split, use_splitting);

        [sorted_locs, sorted_vals, sorted_widths, pairs, error_coeff] = ...
            partner_peak(locs, centers, tolerance, signal, threshold, freq, ...
                vals, widths, proms, method, false, raw);

        centers = get_centers(sorted_locs, freq, split, use_splitting);
        baseline = find_baseline(y, x, data, freq);
        seed = peaks_to_seed(sorted_locs, sorted_vals, sorted_widths, pairs, ...
            baseline, error_coeff + error, seed_method, centers);
        % parameters_dataframe(y, x, :) = seed;
end

function centers = get_centers(locs, freq, split, use_splitting)
    if use_splitting
        centers = zeros(size(locs));
        maxX = max(freq); % If the ideal is past the measured frequency
        for i=1:numel(locs)
            centers(i) = locs(i) + 2 * split;
            if locs(i) + split > maxX
                centers(i) = locs(i) - 2 * split;
            end
        end
    else
        centers = 2 * split - locs; % (split - locs) + split
    end
end
