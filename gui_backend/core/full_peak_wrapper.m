function [pdf, lock_laplacian] = full_peak_wrapper(all_param_data, data, freq, params, lock_laplacian, ref_data, num_params)
    if nargin < 2 && ~isempty(all_param_data)
        data = all_param_data.gWide.signal;
        freq = all_param_data.gWide.SweepParam;
        params = all_param_data.params;
        lock_laplacian = 0;
        ref_data = all_param_data.reflection_df;
        num_params = numel(all_param_data.param_labels) - 1;
    end
    % mpiprofile on;
    f = waitbar(0, "Starting");
    xsize = size(data, 2);
    ysize = size(data, 1);
    pdf = zeros(ysize, xsize, num_params);
    locs_array = zeros(xsize, ysize, 12);
    waitbar(0.33, f, "Finding Peaks");
    parfor x = 1:xsize
        for y = 1:ysize
            shifting = ref_data(y, x);
            [pdf(y, x, :), new_peaks_info, peaks_info] = peak_find_function(x, y, data, freq, params, shifting);
            locs = union(new_peaks_info.locs, peaks_info.locs);
            locs_mat = zeros(12, 1);
            locs_mat(1:numel(locs)) = locs;
            locs_array(y, x, :) = locs_mat;
        end
    end
    % locs_array(:, :, all(locs_array == 0)) = []; % removes column if the entire column is zero
    % save('data.mat', "locs_array")
    waitbar(0.66, f, "Smoothing");
    pdf = smoothing(pdf, params.smoothing_method);
    if params.lock_flag
        [pdf, ~] = laplacian_smoothing(pdf, lock_laplacian, false);
    else
        [pdf, lock] = laplacian_smoothing(pdf, 0, false);
        lock_laplacian = lock;
    end
    if params.smoothing_method == 2
        % Repairs should use the old value if we're locked and
        % make new ones if we're not
        if params.lock_flag
            lock = lock_laplacian;
        else
            lock = 0;
        end
        [pdf, thresh] = laplacian_repair(lock, pdf, ...
            data, freq, ref_data, params, false);
        lock_laplacian = thresh;
    end
    waitbar(1, f, "Completed");
    close(f);
end

