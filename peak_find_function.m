function [seed, threshold, pairs, new_locs, new_vals, centers, locs, vals] = ...
    peak_find_function(use_splitting, split, tolerance, x, y, data, freq, num_max_peaks, peak_perc_threshold, ...
        diff_peak_distance, smooth_span, smooth_degree, method, seed_method, normalize)
        if nargin < 15
            normalize = false;
        end
        if nargin < 8
            num_max_peaks = 2;
            peak_perc_threshold = 5;
            diff_peak_distance = 0.002;
            smooth_span = 41;
            smooth_degree = 7;
            method = 1;
            seed_method = 1;
        end

        raw = squeeze(data(x, y, :));
        if normalize
            raw=raw-mean(raw);
            raw=raw/std(raw);
        end
        signal = smooth(raw, smooth_span, 'sgolay', smooth_degree);
        [vals, locs, widths, proms, threshold, error] = find_peaks_at_point(signal, freq, ...
            false, num_max_peaks, peak_perc_threshold, diff_peak_distance);
 
        % Sort peaks from highest to lowest
        [vals, sorted_indices] = sort(vals);
        sorted_indices = sorted_indices';
        locs = locs(sorted_indices); % Rearrange accordingly
        widths = widths(sorted_indices);
        proms = proms(sorted_indices);
        centers = get_centers(locs, freq, split, use_splitting);

        [new_locs, new_vals, new_widths, pairs, error_coeff] = ...
            partner_peak(locs, centers, tolerance, signal, threshold, freq, ...
                vals, widths, proms, method, false, raw);

        % Get new centers with new locations
        centers = get_centers(new_locs, freq, split, use_splitting);
        baseline = find_baseline(y, x, data, freq);
        seed = peaks_to_seed(new_locs, new_vals, new_widths, pairs, ...
            baseline, error_coeff + error, seed_method, centers);
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
