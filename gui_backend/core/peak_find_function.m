function [seed, new_peaks_info, peaks_info] = ...
    peak_find_function(x, y, data, freq, params_struct, shift)
        if nargin > 5
            params_struct.shift = shift;
        end

        raw = squeeze(data(x, y, :));
        if params_struct.normalize
            raw=raw-mean(raw);
            raw=raw/std(raw);
        end
        signal = smooth(raw, params_struct.smooth_span, ...
                        'sgolay', params_struct.smooth_degree);
        
        peaks_info = find_peaks_at_point(signal, freq, false, params_struct);
        peaks_info.centers = get_centers(peaks_info, params_struct, freq);
        remove_idx = filter_peaks(peaks_info, params_struct);
        peaks_info.vals(remove_idx) = [];
        peaks_info.locs(remove_idx) = [];
        peaks_info.widths(remove_idx) = [];
        peaks_info.proms(remove_idx) = [];
        peaks_info.centers(remove_idx) = [];

        gen_peaks_info = partner_peak(signal, raw, freq, peaks_info, params_struct, false);
        % Regenerate centers for generated peaks
        new_peaks_info = gen_peaks_info;
        new_peaks_info.centers = get_centers(new_peaks_info, params_struct, freq);
        duplicates = new_peaks_info.pairs(new_peaks_info.pairs(:, 1) == new_peaks_info.pairs(:, 2));
        new_peaks_info.pairs(duplicates, :) = [];
        new_peaks_info.centers(duplicates) = [];
        new_peaks_info.error = new_peaks_info.error + peaks_info.error;
        new_peaks_info.baseline = find_baseline(y, x, data, freq);
        
        seed = peaks_to_seed(new_peaks_info, params_struct);
end

function centers = get_centers(peaks_info, params_struct, freq)
    locs = peaks_info.locs;
    partnering_method = params_struct.partnering_method;
    shift = params_struct.shift;
    centers = zeros(size(locs));
    if partnering_method == 1
        centers = 2 * shift - locs; % (shift - locs) + shift
    elseif partnering_method == 2
        maxX = max(freq); % If the ideal is past the measured frequency
        for i=1:numel(locs)
            centers(i) = locs(i) + 2 * shift;
            if locs(i) + shift > maxX
                centers(i) = locs(i) - 2 * shift;
            end
        end
    elseif partnering_method == 3
        centers = flip(locs);
    end
end
