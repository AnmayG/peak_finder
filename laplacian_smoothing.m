% Takes in a parameters_dataframe matrix including the E map and tries
% Then computes the Laplacian and goes through by number of tries
% Idea is that the harder we had to try, the more likely we got it wrong
% Thus, we need to go through all the values that we tried and see how bad
% the Laplacian is, then fix it pixel-by-pixel

% Idea for this stage is that the Laplacian acts as an error metric telling
% us how wrong the guess is according to continuity, so by minimizing the
% Laplacian we can get some

function parameters_dataframe = laplacian_smoothing(input_params, signal, freq, split, tolerance, num_max_peaks, ...
                            peak_perc_threshold, diff_peak_distance, smooth_span, smooth_degree, display)
    if nargin < 4
        display = false;
    end
    parameters_dataframe = input_params;
    lpeak = squeeze(input_params(:, :, 15)) / 1e6;
    rpeak = squeeze(input_params(:, :, 16)) / 1e6;
    tries = squeeze(input_params(:, :, 14));

    e111 = squeeze(parameters_dataframe(:, :, 2));
    % Express positions in megahertz for ease of display
    % Also because it makes understanding the limits easier 
    % (a 15 megahertz jump is probably unphysical...)
    e111 = e111 ./ 1e6;
    laplacian = del2(e111);
    output = zeros(size(laplacian, 1), size(laplacian, 2));
    perc_thresholds = [10, 34, 68, 90];
    lapl_thresholds = prctile(laplacian, perc_thresholds, 'all');
    
    xsize = size(laplacian, 2);
    ysize = size(laplacian, 1);
    for x=1:xsize
        for y=1:ysize
            % if tries(x, y) ~= 0
                lapl = laplacian(x, y);
                sig = squeeze(signal(x, y, :));
                if lapl < lapl_thresholds(1) % Under 5%
                    % Harsh fix, splittings are too close
                    % Complete failure to fit, mark it red?
                    % if trial == 1
                    output(x, y) = -2;
                    % else
                    %     [vals, locs, widths, proms, threshold, signal2, raw] = ...
                    %         find_peaks_at_point(y, x, signal, freq, false, num_max_peaks, ...
                    %         peak_perc_threshold, diff_peak_distance, smooth_span, smooth_degree);
                    %     [sorted_locs, sorted_vals, sorted_widths, pairs, error_coeff] = ...
                    %         partner_peak(locs, split, tolerance, signal2, threshold, freq, vals, widths, proms, 3, false, raw);
                    %     baseline = find_baseline(y, x, data, freq);
                    %     seed = peaks_to_seed(sorted_locs, sorted_vals, sorted_widths, pairs, baseline, error_coeff);
                    %     parameters_dataframe(x, y, :) = seed;
                    % end
                elseif lapl < lapl_thresholds(2) % 5-34%
                    % Medium increase
                    output(x, y) = 0;
                elseif lapl < lapl_thresholds(3) % 34-68%
                    % Good enough/slight fix
                    output(x, y) = 0;
                elseif lapl < lapl_thresholds(4) % 68-95%
                    % Medium decrease
                    output(x, y) = 0;
                else % Over 5%
                    % Harsh fix, splittings are too far
                    % Complete failure to fit, mark it red?
                    output(x, y) = 2;
                end
            % end
        end
    end
    whos parameters_dataframe
    parameters_dataframe(:, :, end + 1) = output; 
    whos parameters_dataframe

    if display
        figure;
        subplot(2, 3, 1);
        hold on;
        e111(e111 == 0) = NaN;
        set(gca, 'color', [1 0 0]);
        imagesc(parameters_dataframe(:, :, 2), 'AlphaData', isnan(e111)); % Display the slice
        colorbar; % Add a colorbar to each subplot
        hold on;
        colorbar;
        clim('auto');
        axis tight;
        axis equal;
        title("E_{111}");
        hold off;

        subplot(2, 3, 4);
        hold on;
        e111(e111 == 0) = NaN;
        set(gca, 'color', [1 0 0]);
        imagesc(parameters_dataframe(:, :, 15), 'AlphaData', isnan(e111)); % Display the slice
        colorbar; % Add a colorbar to each subplot
        clim('auto');
        axis tight;
        axis equal;
        title("Left Peak");
        hold off;
        
        subplot(2, 3, 5);
        hold on;
        e111(e111 == 0) = NaN;
        set(gca, 'color', [1 0 0]);
        imagesc(parameters_dataframe(:, :, 16), 'AlphaData', isnan(e111)); % Display the slice
        colorbar; % Add a colorbar to each subplot
        clim('auto');
        axis tight;
        axis equal;
        title("Right Peak");
        hold off;
        
        subplot(2, 3, 2);
        hold on;
        imagesc(tries);
        colorbar;
        clim('auto');
        axis tight;
        axis equal;
        title("Tries");
        hold off;

        subplot(2, 3, 3);
        hold on;
        imagesc(laplacian);
        colorbar;
        clim('auto');
        axis tight;
        axis equal;
        title("Laplacian");
        hold off;
    end
end