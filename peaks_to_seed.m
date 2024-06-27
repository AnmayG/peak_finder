function parameters = peaks_to_seed(locs, vals, widths, pairs, ...
    baseline, error_coeff, method, guess)
    if nargin < 7
        method = 1;
    end

    % Format into parameters
    parameters = zeros(16, 1);

    if ~isempty(pairs)
        peaks_111_ind = decide_peaks(locs, vals, pairs, guess, method);
        peaks_111 = pairs(peaks_111_ind, :);
        peaks_111 = peaks_111(1, :);
        % Locations should be sorted
        if locs(peaks_111(1)) > locs(peaks_111(2))
            peak1 = peaks_111(1);
            peaks_111(1) = peaks_111(2);
            peaks_111(2) = peak1;
        end
        
        try
            peaks_111_locs = [locs(peaks_111(1)), locs(peaks_111(2))];
            peaks_111_vals = [vals(peaks_111(1)), vals(peaks_111(2))];
            peaks_111_fwhm = [widths(peaks_111(1)), widths(peaks_111(2))];
            peaks_111_cont = [(peaks_111_vals(1) - baseline) / baseline, (peaks_111_vals(2) - baseline) / baseline];
            peaks_111_stds = peaks_111_fwhm ./ (2 * sqrt(2 * log(2))); % assuming Gaussian
        catch
            % :(
        end
        
        parameters(1) = (peaks_111_locs(1) + peaks_111_locs(2)) / 2; % Average split
        parameters(2) = abs(peaks_111_locs(2) - peaks_111_locs(1)) / 2; % Deviation
        parameters(3) = (peaks_111_stds(1) + peaks_111_stds(2)) / 2; % Average width
        parameters(4) = abs(peaks_111_stds(1) - peaks_111_stds(2)) / 2; % Deviation
        parameters(5) = (peaks_111_cont(1) + peaks_111_cont(2)) / 2; % Average height
        parameters(6) = (peaks_111_cont(2) - peaks_111_cont(1)) / 2; % Deviation
        
        pairs(peaks_111_ind, :) = [];
        guess(peaks_111_ind) = [];
        if ~isempty(pairs)
            peak_height_sums = decide_peaks(locs, vals, pairs, guess, method);
            [~, peaks_non111_ind] = min(peak_height_sums);
            peaks_non111 = pairs(peaks_non111_ind, :);
            peaks_non111_locs = [locs(peaks_non111(1)), locs(peaks_non111(2))];
            peaks_non111_vals = [vals(peaks_non111(1)), vals(peaks_non111(2))];
            peaks_non111_fwhm = [widths(peaks_non111(1)), widths(peaks_non111(2))];
            peaks_non111_cont = [(peaks_non111_vals(1) - baseline) / baseline, (peaks_non111_vals(2) - baseline) / baseline];
            peaks_non111_stds = peaks_non111_fwhm ./ (2 * sqrt(2 * log(2))); % assuming Gaussian
    
            parameters(7) = (peaks_non111_locs(1) + peaks_non111_locs(2)) / 2; % Average split
            parameters(8) = abs(peaks_non111_locs(2) - peaks_non111_locs(1)) / 2; % Deviation
            parameters(9) = (peaks_non111_stds(1) + peaks_non111_stds(2)) / 2; % Average width
            parameters(10) = abs(peaks_non111_stds(1) - peaks_non111_stds(2)) / 2; % Deviation
            parameters(11) = peaks_non111_cont(1) + peaks_non111_cont(2); % Average height
            parameters(12) = (peaks_non111_cont(2) - peaks_non111_cont(1)) / 2; % Deviation
        end

        parameters(13) = baseline;

        parameters(14) = error_coeff;
        parameters(15) = min(peaks_111_locs);
        parameters(16) = max(peaks_111_locs);
    else
        % Can use global parameters here if we want, this is a failure case
        % parameters(1) = 2.7e9;
        % parameters(7) = 1;
    end
end

function peaks_111_ind = decide_peaks(locs, vals, pairs, guess, method)
    if method == 1 % Tallest single
        peak_height_sums = min(vals(pairs(:, 1)), vals(pairs(:, 2)));
        [~, peaks_111_ind] = min(peak_height_sums);
    elseif method == 2 % Tallest avg
        peak_height_sums = vals(pairs(:, 1)) + vals(pairs(:, 2)); % tallest average
        [~, peaks_111_ind] = min(peak_height_sums);
    elseif method == 3 % Tallest out
        [~, leftmost] = min(min(locs(pairs(:, 1)), locs(pairs(:, 2))));
        [~, rightmost] = max(max(locs(pairs(:, 1)), locs(pairs(:, 2))));
        outermost = [leftmost; rightmost];
        vals(setdiff(1:end, outermost)) = Inf;
        peak_height_sums = min(vals(pairs(:, 1)), vals(pairs(:, 2)));
        [~, peaks_111_ind] = min(peak_height_sums);
    elseif method == 4 % Leftmost
        peak_height_sums = min(locs(pairs(:, 1)), locs(pairs(:, 2)));
        [~, peaks_111_ind] = min(peak_height_sums);
    elseif method == 5 % Rightmost
        peak_height_sums = max(locs(pairs(:, 1)), locs(pairs(:, 2)));
        [~, peaks_111_ind] = max(peak_height_sums);
    elseif method == 6 % Innermost - smallest splitting
        peak_height_sums = abs(locs(pairs(:, 1)) - locs(pairs(:, 2)));
        [~, peaks_111_ind] = min(peak_height_sums);
    elseif method == 7 % closest to guess
        peak_height_sums = (locs(pairs(:, 1)) + locs(pairs(:, 2))) / 2; % Get the shift
        peak_height_sums = abs(peak_height_sums - guess); % Distance from guessed centers
        [~, peaks_111_ind] = min(peak_height_sums);
    else
        peak_height_sums = min(vals(pairs(:, 1)), vals(pairs(:, 2))); % tallest total
        [~, peaks_111_ind] = min(peak_height_sums);
    end
end
