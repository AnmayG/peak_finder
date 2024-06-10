function parameters_dataframe = smoothing(parameters_dataframe, method)
    if method == 1
        return;
    elseif method == 2
        % Laplacian fix
    elseif method == 3
        % Average
        K = (1/16)*ones(4);
        for index = 1:size(parameters_dataframe, 3)
            parameters_dataframe(:, :, index) = conv2(parameters_dataframe(:, :, index), K, 'same');
        end
    elseif method == 4
        % Gaussian
        for index = 1:size(parameters_dataframe, 3)
            parameters_dataframe(:, :, index) = imgaussfilt(parameters_dataframe(:, :, index), 2);
        end
    end
end

