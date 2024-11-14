function remove_idx = filter_peaks(peaks_info, params_struct, old_loc)
    if nargin < 3 && ismember(params_struct.zoning_method, [4, 5, 7])
        % No old location provided, not zoning
        old_loc = -1;
        % params_struct.zoning_method = 1;
    end
    if(isfield(params_struct, "zones_x"))
        zones_x = params_struct.zones_x;
        zones_y = params_struct.zones_y;
        zoning_method = params_struct.zoning_method;
    else
        zones_x = [0, Inf];
        zones_y = [0, Inf];
        zoning_method = 1;
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
            if(old_loc == -1)
                return;
            end
            for i=1:length(locs)
                shift = (old_loc + locs(i)) / 2;
                if ~is_valid_peak(shift, zones_x)
                    remove_idx(i) = true;
                end
            end
        case 5
            % Split
            if(old_loc == -1)
                return;
            end
            for i=1:length(locs)
                split = (old_loc - locs(i)) / 2;
                if ~is_valid_peak(split, zones_y)
                    remove_idx(i) = true;
                end
            end
            disp(remove_idx)
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
            if(old_loc == -1)
                return;
            end
            freq_zones = zones_x;
            for i=1:2:numel(zones_x) - 2
                disp(i)
                disp(zones_y)
                freq_zones(i) = zones_x(i) - zones_y(i);
                freq_zones(i+1) = zones_x(i+1) + zones_y(i+1);
            end
            for i=1:length(locs)
                if ~is_valid_peak(locs(i), freq_zones)
                    remove_idx(i) = true;
                end
            end
            if old_loc ~= -1
                for i=1:length(locs)
                    shift = (old_loc + locs(i)) / 2;
                    if ~is_valid_peak(shift, zones_x)
                        remove_idx(i) = true;
                    end
                end
                for i=1:length(locs)
                    split = (old_loc - locs(i)) / 2;
                    if ~is_valid_peak(split, zones_y)
                        remove_idx(i) = true;
                    end
                end
            end
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
