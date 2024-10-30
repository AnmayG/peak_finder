function remove_idx = filter_peaks(peaks_info, params_struct)
    if(isfield(params_struct, "zones_x"))
        zones_x = params_struct.zones_x;
        zones_y = params_struct.zones_y;
        zoning_method = params_struct.zoning_method;
    else
        zones_x = [0, Inf];
        zones_y = [0, Inf];
        zoning_method = 0;
    end
    % Zoning filters
    % If we're doing absolute zoning, peaks outside the allowed values
    % should be immediately removed.
    locs = peaks_info.locs;
    vals = peaks_info.vals;
    remove_idx = false(length(locs));
    switch zoning_method
        case 2
            % Frequency
            for i=1:length(locs)
                if ~is_valid_peak(locs(i), zones_x)
                 remove_idx(i) = true;
                end
            end
        case 3
            % Values
            for i=1:length(vals)
                if ~is_valid_peak(vals(i), zones_y)
                 remove_idx(i) = true;
                end
            end
        case 4
            % Shift
        case 5
            % Split
        case 6
            % Frequency + Vals
            for i=1:length(locs)
                if ~is_valid_peak(locs(i), zones_x)
                 remove_idx(i) = true;
                end
            end
            for i=1:length(vals)
                if ~is_valid_peak(vals(i), zones_y)
                 remove_idx(i) = true;
                end
            end
        case 7
            % Split + Shift
    end
end

function result = is_valid_peak(coord, valid_locs)
    result = false;
    num_intervals = length(valid_locs) - 1;
    for i = 1:2:num_intervals
        if coord >= valid_locs(i) && coord <= valid_locs(i+1)
            result = true;
            return;
        end
    end
end
