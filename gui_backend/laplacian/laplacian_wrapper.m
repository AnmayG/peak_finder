function [pdf, lock] = laplacian_wrapper(params, pdf, peak_info_df, lock_laplacian)
    % make new ones if we're not
    if params.lock_flag
        lock = lock_laplacian;
    else
        lock = 0;
    end
    for i=1:500
        % laplacian_repair does not change peak_info
        [pdf2, thresh] = laplacian_repair(lock, pdf, peak_info_df, 200);
        lock_laplacian = thresh;
        lock = lock_laplacian;
        if ~isequal(squeeze(pdf2(:, :, 17)), squeeze(pdf(:, :, 17)))
            pdf = pdf2;
        else
            break
        end
    end
end

