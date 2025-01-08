function peaks_info = ...
    find_peaks_at_point(signal, freq, display, params_struct)
        if nargin < 3 % Don't need to include display if not necessary
            display = false;
            params_struct = struct();
            params_struct.num_max_peaks = 2;
            % Previously 32%
            params_struct.peak_perc_threshold = 5;
            params_struct.diff_peak_distance = 0.02e9;
        end
        error = 0;
        x_bounds = (freq > freq(1)) & (freq < freq(end));
        peak_threshold = mean(signal) - params_struct.peak_perc_threshold * std(signal);
        norm_signal = max(-signal + peak_threshold, 0);
        adj_signal = norm_signal(x_bounds);
        adj_freq = freq(x_bounds);
        
        % Zone splicing
        % Note that currently, zones are naively spliced together
        % If a frequency is out of range, I don't use it
        mask1 = false(size(adj_freq));
        if params_struct.zoning_method == 2 || params_struct.zoning_method == 6
            zones = reshape(params_struct.zones_x, 2, []).';
            comparisons = (adj_freq >= zones(:,1)) & (adj_freq <= zones(:,2));
            mask1 = ~any(comparisons, 1);
        end
        mask2 = false(size(adj_freq));
        if params_struct.zoning_method == 3 || params_struct.zoning_method == 6
            zones = reshape(params_struct.zones_y, 2, []).';
            comparisons = (signal(x_bounds) >= zones(:,1)) & (signal(x_bounds) <= zones(:,2));
            mask2 = ~any(comparisons, 2)';
        end
        if params_struct.zoning_method == 6
            mask = mask1 & mask2;
        elseif params_struct.zoning_method == 7
            freq_zones = params_struct.zones_x;
            for i=1:2:length(params_struct.zones_x)
                if i > length(params_struct.zones_x) || i > length(params_struct.zones_y)
                    break;
                end
                freq_zones(i) = params_struct.zones_x(i) - params_struct.zones_y(i);
                freq_zones(i+1) = params_struct.zones_x(i+1) + params_struct.zones_y(i+1);
            end
            zones = reshape(freq_zones, 2, []).';
            comparisons = (adj_freq >= zones(:,1)) & (adj_freq <= zones(:,2));
            mask = ~any(comparisons, 1);
        else
            mask = mask1 | mask2;
        end
        adj_signal(mask) = [];
        adj_freq(mask) = [];

        [vals, locs, widths, proms] = findpeaks_custom(adj_signal, adj_freq, ...
                            'SortStr','descend',...
                            'WidthReference','halfheight', ...
                            'NPeaks', params_struct.num_max_peaks, ...,
                            'MinPeakDistance', params_struct.diff_peak_distance);
        vals = -vals + peak_threshold;
        % Somehow these get transposed idk why
        vals = vals';
        proms = proms';
        % Sort peaks from left to right
        [~, sorted_indices] = sort(locs);
        sorted_indices = sorted_indices';
        locs = locs(sorted_indices)';
        vals = vals(sorted_indices)'; % Rearrange accordingly
        widths = widths(sorted_indices)';
        proms = proms(sorted_indices)';
        
        peaks_info = struct();
        peaks_info.vals = vals;
        peaks_info.locs = locs;
        peaks_info.widths = widths;
        peaks_info.proms = proms;
        peaks_info.threshold = peak_threshold;
        peaks_info.error = error;

        if params_struct.refine
            [~, indices] = arrayfun(@(x) min(abs(adj_freq - x)), locs); %extract peak indicies

            [~, peaks_info.locs]=refinepeaks(adj_signal, indices, adj_freq); %toss the updated vals since there is a scaling issue and they're mostly aesthetic/hidden by the size of the peak marker
        end


        if display
            figure
            hold on
            yline(peak_threshold, 'Color', 'r')
            plot(freq, signal, 'Color', 'b')
            plot(locs, vals, 'rv', 'MarkerFaceColor', 'g');
            title("Point at " + x + ", " + y)
            disp(x + " " + y + " completed");
        end
end
