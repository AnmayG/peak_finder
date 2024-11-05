function [pdf, lock_laplacian, peak_info_df] = full_peak_wrapper(all_param_data, data, freq, params, lock_laplacian, ref_data, num_params)
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
    
    max_num_peaks = 16;
    locs_array = zeros(xsize, ysize, max_num_peaks);
    vals_array = zeros(xsize, ysize, max_num_peaks);
    widths_array = zeros(xsize, ysize, max_num_peaks);
    proms_array = zeros(xsize, ysize, max_num_peaks);
    baseline_array = zeros(xsize, ysize);
    mask = params.roi_mask;

    waitbar(0.33, f, "Finding Peaks");
    parfor x = 1:xsize
        for y = 1:ysize
            if ~mask(y,x)
                pdf(y, x, :) = zeros(num_params, 1);
                locs_array(y, x, :) = zeros(1, max_num_peaks);
                vals_array(y, x, :) = zeros(1, max_num_peaks);
                widths_array(y, x, :) = zeros(1, max_num_peaks);
                proms_array(y, x, :) = zeros(1, max_num_peaks);
                baseline_array(y, x) = 0;
                disp([y, x]);
                continue;
            end
            shifting = ref_data(y, x);
            [pdf(y, x, :), new_peaks_info, peaks_info] = peak_find_function(x, y, data, freq, params, shifting);
            
            % To an end user, unpartnered/generated peaks are still useful
            % info, so when we return locs/vals/widths/proms it should be
            % union based
            [locs, ia, ib] = union(new_peaks_info.locs, peaks_info.locs);
            vals = [new_peaks_info.vals(ia)', peaks_info.vals(ib)'];
            widths = [new_peaks_info.widths(ia)', peaks_info.widths(ib)'];
            proms = [new_peaks_info.proms(ia)', peaks_info.proms(ib)'];
            
            locs_mat = zeros(max_num_peaks, 1);
            locs_mat(1:numel(locs)) = locs;
            locs_array(y, x, :) = locs_mat;
            vals_mat = zeros(max_num_peaks, 1);
            vals_mat(1:numel(vals)) = vals;
            vals_array(y, x, :) = vals_mat;
            widths_mat = zeros(max_num_peaks, 1);
            widths_mat(1:numel(widths)) = widths;
            widths_array(y, x, :) = widths_mat;
            proms_mat = zeros(max_num_peaks, 1);
            proms_mat(1:numel(proms)) = proms;
            proms_array(y, x, :) = proms_mat;
            baseline_array(y, x) = new_peaks_info.baseline;
        end
    end
    zero_slices = squeeze(all(all(locs_array == 0, 1), 2));
    locs_array(:, :, zero_slices) = []; % removes column if the entire column is zero
    vals_array(:, :, zero_slices) = [];
    widths_array(:, :, zero_slices) = [];
    proms_array(:, :, zero_slices) = [];
    peak_info_df = struct("locs", locs_array, "vals", vals_array, ...
        "widths", widths_array, "proms", proms_array, "baselines", baseline_array);
    waitbar(0.66, f, "Smoothing");
    pdf = smoothing(pdf, params.smoothing_method);
    if params.lock_flag
        [pdf, ~] = laplacian_calculate(pdf, lock_laplacian, false);
    else
        [pdf, lock] = laplacian_calculate(pdf, 0, false);
        lock_laplacian = lock;
    end
    if params.smoothing_method == 2
        [pdf, lock_laplacian] = laplacian_wrapper(params, pdf, peak_info_df);
    end
    waitbar(1, f, "Completed");
    close(f);
end

