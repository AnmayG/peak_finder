load('327_s1_20p1gpa_126K_mwn48_z6A.mat');

x = 1;
y = 1;
x_grid_width = 60;
y_grid_width = 60;
xsize = 180;
ysize = 180;
split = 3.1e9;
tolerance = 0.1e9;
for x = randi([0, x_grid_width - 1], 1, xsize / x_grid_width) + (1:x_grid_width:xsize)
    for y = randi([0, y_grid_width - 1], 1, ysize / y_grid_width) + (1:y_grid_width:ysize)    
        [vals, locs, widths, ~, threshold, signal, raw] = find_peaks_at_point(x, y, gWide);
        [sorted_locs, sorted_vals, sorted_widths, pairs] = partner_peak(locs, split, tolerance, signal, threshold, gWide, vals, widths, true, raw);
        baseline = find_baseline(x, y, gWide);
        seed = peaks_to_seed(sorted_locs, sorted_vals, sorted_widths, pairs, baseline, split, tolerance);
    end
end
