% Takes in a parameters_dataframe matrix including the E map and tries
% Then computes the Laplacian and goes through by number of tries
% Idea is that the harder we had to try, the more likely we got it wrong
% Thus, we need to go through all the values that we tried and see how bad
% the Laplacian is, then fix it pixel-by-pixel

% Idea for this stage is that the Laplacian acts as an error metric telling
% us how wrong the guess is according to continuity, so by minimizing the
% Laplacian we can get some

function [parameters_dataframe, lapl_thresholds] = laplacian_calculate(input_params, lock, display)
    parameters_dataframe = input_params;

    e111 = squeeze(parameters_dataframe(:, :, 2));
    e111 = e111 ./ 1e6;
    laplacian = del2(e111);
    output = zeros(size(laplacian, 1), size(laplacian, 2));
    perc_thresholds = [12, 34, 68, 88];
    if lock ~= 0 % Lock should be 0 if it doesn't exist
        lapl_thresholds = lock;
    else
        lapl_thresholds = prctile(laplacian, perc_thresholds, 'all');
    end

    xsize = size(laplacian, 2);
    ysize = size(laplacian, 1);
    for x=1:xsize
        for y=1:ysize
            lapl = laplacian(x, y);
            if lapl < lapl_thresholds(1) % Under 5%
                % Harsh fix, splittings are too close
                % Complete failure to fit, mark it red
                output(x, y) = 1;
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
                % Complete failure to fit, mark it red
                output(x, y) = 1;
            end
        end
    end

    % Fill in the gaps
    se = strel("disk", 1);
    output = imclose(output, se);
    output = imfill(output, 8);
    parameters_dataframe(:, :, 17) = output;
end