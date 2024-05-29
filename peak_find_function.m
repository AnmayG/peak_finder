function seed = peak_find_function(split, tolerance, x, y, gWide)
        [vals, locs, widths, ~, threshold, signal, raw] = find_peaks_at_point(x, y, gWide, false, 2, 5, 1e8, 41, 6);
        [sorted_locs, sorted_vals, sorted_widths, pairs] = partner_peak(locs, split, tolerance, signal, threshold, gWide, vals, widths, false, raw);
        baseline = find_baseline(x, y, gWide);
        seed = peaks_to_seed(sorted_locs, sorted_vals, sorted_widths, pairs, baseline, split, tolerance);
end
