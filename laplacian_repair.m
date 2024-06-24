% Takes in a parameters_dataframe matrix including the E map and tries
% Then computes the Laplacian and goes through by number of tries
% Idea is that the harder we had to try, the more likely we got it wrong
% Thus, we need to go through all the values that we tried and see how bad
% the Laplacian is, then fix it pixel-by-pixel

% Idea for this stage is that the Laplacian acts as an error metric telling
% us how wrong the guess is according to continuity, so by minimizing the
% Laplacian we can get some

function [parameters_dataframe, thresh] = laplacian_repair(use_splitting, lock, input_params, data, freq, splits, tolerance, num_max_peaks, ...
                            peak_perc_threshold, diff_peak_distance, smooth_span, smooth_degree, seed_method, display)
    parameters_dataframe = input_params;
    % lpeak = squeeze(input_params(:, :, 15)) / 1e6;
    % rpeak = squeeze(input_params(:, :, 16)) / 1e6;
    tries = squeeze(input_params(:, :, 14));

    e111 = squeeze(parameters_dataframe(:, :, 2));
    laplacian = squeeze(parameters_dataframe(:, :, 17));

    xsize = size(laplacian, 2);
    ysize = size(laplacian, 1);
    if isscalar(splits)
        splits = splits * ones(xsize, ysize);
    end

    parfor x=1:xsize
        for y=1:ysize
            if laplacian(y, x) == 1 % If flagged by laplacian
                split = splits(y, x);
                parameters_dataframe(y, x, :) = [peak_find_function(use_splitting, split, ...
                    tolerance, x, y, data, freq, num_max_peaks, peak_perc_threshold, ...
                    diff_peak_distance, smooth_span, smooth_degree, 3, seed_method); 0]; % Reseed with closest
            end
        end
    end
    [parameters_dataframe, thresh] = laplacian_smoothing(parameters_dataframe, lock, false);
    
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