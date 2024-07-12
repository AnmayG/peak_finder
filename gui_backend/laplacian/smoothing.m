function parameters_dataframe = smoothing(parameters_dataframe, method)
    if method == 1
        return;
    elseif method == 2
        % Laplacian fix
    elseif method == 3
        % Average
        sizes = 4;
        K = (1/sizes^2)*ones(sizes);
        for index = 1:size(parameters_dataframe, 3)
            parameters_dataframe(:, :, index) = conv2(parameters_dataframe(:, :, index), K, 'same');
        end
    elseif method == 4
        % Gaussian
        for index = 1:size(parameters_dataframe, 3)
            df = parameters_dataframe(:, :, index);
            df(isnan(df)) = 0;
            parameters_dataframe(:, :, index) = imgaussfilt(df, 3);
        end
    elseif method == 5
        % Spatially variant Laplacian blur
        for index = 1:size(parameters_dataframe, 3)
            df = parameters_dataframe(:, :, index);
            df(isnan(df)) = 0;
            
            parameters_dataframe(:, :, index) = imgaussfilt(df, 3);
        end
    end
end

