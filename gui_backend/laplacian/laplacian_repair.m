% Takes in a parameters_dataframe matrix including the E map and tries
% Then computes the Laplacian and goes through by number of tries
% Idea is that the harder we had to try, the more likely we got it wrong
% Thus, we need to go through all the values that we tried and see how bad
% the Laplacian is, then fix it pixel-by-pixel
function [parameters_dataframe, thresh] = laplacian_repair(lock, input_params, peaks_info_array, radius)
    parameters_dataframe = input_params;
    laplacian = squeeze(parameters_dataframe(:, :, 17));

    xsize = size(laplacian, 2);
    ysize = size(laplacian, 1);
    for x=1:xsize
        for y=1:ysize
            if laplacian(y, x) == 1 % If flagged by laplacian
                pixel_params = minimize_laplacian(y, x, peaks_info_array, parameters_dataframe, radius - 1);
                parameters_dataframe(y, x, :) = pixel_params;
            end
        end
    end
    [parameters_dataframe, thresh] = laplacian_calculate(parameters_dataframe, lock, false);
end

function [pixel_params] = minimize_laplacian(y, x, peaks_info_array, input_params, radius)
    locs_array = nonzeros(squeeze(peaks_info_array.locs(y, x, :)));
    vals_array = nonzeros(squeeze(peaks_info_array.vals(y, x, :)));
    widths_array = nonzeros(squeeze(peaks_info_array.widths(y, x, :)));
    baseline = squeeze(peaks_info_array.baselines(y, x));
    if numel(locs_array) < 2
        pixel_params = input_params(y, x, :);
        return
    end
    
    min_laplacian = inf;
    mini = 1;
    minj = 1;
    peaks_info = struct("locs", [], "vals", [], "widths", [], "baseline", baseline);
    neighbors = del_at_point(y, x, squeeze(input_params(:, :, 2)), radius);
    for i=1:size(locs_array, 1)
        for j=i:size(locs_array, 1)
            if i == j
                continue
            end
            curr_val = abs(locs_array(j) - locs_array(i)) / 2;
            lapl = abs(neighbors - curr_val);
            if lapl < min_laplacian
                peaks_info.locs = locs_array([i, j]);
                peaks_info.vals = vals_array([i, j]);
                peaks_info.widths = widths_array([i, j]);
                mini = i;
                minj = j;
                min_laplacian = lapl;
            end
        end
    end
    % if y == 121 && x == 63
    %     disp(locs_array')
    %     disp([mini, minj])
    %     disp(min_laplacian)
    % end

    if numel(locs_array) > 2
        locs_array([mini, minj]) = [];
        vals_array([mini, minj]) = [];
        widths_array([mini, minj]) = [];
        min_laplacian = inf;
        non111_peaks_info = struct("locs", [], "vals", [], "widths", []);
        neighbors = del_at_point(y, x, squeeze(input_params(:, :, 8)), radius);
        for i=1:size(locs_array, 1)
            for j=1:size(locs_array, 1)
                if i == j
                    continue
                end
                lapl = abs(neighbors - abs(locs_array(j) - locs_array(i)) / 2);
                % shift = abs(neighbors - abs(locs_array(j) + locs_array(i)) / 2);
                if lapl < min_laplacian
                    non111_peaks_info.locs = locs_array([i, j]);
                    non111_peaks_info.vals = vals_array([i, j]);
                    non111_peaks_info.widths = widths_array([i, j]);
                    min_laplacian = lapl;
                end
            end
        end
    else
        non111_peaks_info.locs = [];
        non111_peaks_info.vals = [];
        non111_peaks_info.widths = [];
    end

    pixel_params = [calculate_params(peaks_info, non111_peaks_info); 0];
end

function lapl = del_at_point(y, x, matrix, radius)
    lapl = 0;
    sizes = size(matrix, [1, 2]);
    total = 0;
    for i=max(y-radius, 1):min(y+radius, sizes(1))
        for j=max(x-radius, 1):min(x+radius, sizes(2))
            if i == y && x == j
                continue
            end
            lapl = lapl + matrix(i, j);
            total = total + 1;
        end
    end
    % lapl = lapl / ((radius * 2 + 1)^2 - 1);
    lapl = lapl / total;
end

function parameters = calculate_params(peaks_info, non111_peaks_info)
    parameters = zeros(16, 1);
    baseline = peaks_info.baseline;
    parameters(13) = baseline;
    parameters(14) = 9; % 9 means replaced by Laplacian

    if isempty(peaks_info.vals) || numel(peaks_info.vals) < 2
        return;
    end
    peaks_111_locs = peaks_info.locs;
    peaks_111_vals = peaks_info.vals;
    peaks_111_fwhm = peaks_info.widths;

    peaks_111_cont = [(peaks_111_vals(1) - baseline) / baseline, (peaks_111_vals(2) - baseline) / baseline];
    peaks_111_stds = peaks_111_fwhm ./ (2 * sqrt(2 * log(2))); % assuming Gaussian

    parameters(1) = (peaks_111_locs(1) + peaks_111_locs(2)) / 2; % Average split
    parameters(2) = abs(peaks_111_locs(2) - peaks_111_locs(1)) / 2; % Deviation
    parameters(3) = (peaks_111_stds(1) + peaks_111_stds(2)) / 2; % Average width
    parameters(4) = abs(peaks_111_stds(1) - peaks_111_stds(2)) / 2; % Deviation
    parameters(5) = peaks_111_cont(1) + peaks_111_cont(2); % Average height
    parameters(6) = (peaks_111_cont(2) - peaks_111_cont(1)) / 2; % Deviation
    
    if ~isempty(non111_peaks_info.locs)
        peaks_non111_locs = non111_peaks_info.locs;
        peaks_non111_vals = non111_peaks_info.vals;
        peaks_non111_fwhm = non111_peaks_info.widths;
        peaks_non111_cont = [(peaks_non111_vals(1) - baseline) / baseline, (peaks_non111_vals(2) - baseline) / baseline];
        peaks_non111_stds = peaks_non111_fwhm ./ (2 * sqrt(2 * log(2))); % assuming Gaussian

        parameters(7) = (peaks_non111_locs(1) + peaks_non111_locs(2)) / 2; % Average split
        parameters(8) = abs(peaks_non111_locs(2) - peaks_non111_locs(1)) / 2; % Deviation
        parameters(9) = (peaks_non111_stds(1) + peaks_non111_stds(2)) / 2; % Average width
        parameters(10) = abs(peaks_non111_stds(1) - peaks_non111_stds(2)) / 2; % Deviation
        parameters(11) = peaks_non111_cont(1) + peaks_non111_cont(2); % Average height
        parameters(12) = (peaks_non111_cont(2) - peaks_non111_cont(1)) / 2; % Deviation
    end

    parameters(15) = min(peaks_111_locs);
    parameters(16) = max(peaks_111_locs);
end