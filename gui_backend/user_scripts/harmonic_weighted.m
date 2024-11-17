function function_outputs = harmonic_weighted(function_inputs)
    % Initialize function_outputs as the first input struct
    function_outputs = function_inputs(1);
    parameters_dataframe = function_inputs(1).parameters;
    in = function_inputs(1);
    mask = in.new_settings(1).roi_mask;
    curr_param = function_inputs(1).parameters;

    % Initialize the interpolated data
    pdf = parameters_dataframe; % Copy the original data

    % Only process slice at index = 2
    index = 2;

    % Extract the current slice
    df = parameters_dataframe(:, :, index);

    % Find known points
    [x, y] = meshgrid(1:size(df, 2), 1:size(df, 1));
    known_points = df ~= 0;  % Known points (non-zero values)

    % If no known points, skip this slice
    if any(known_points, 'all')
        % Initialize the harmonic solution with the input slice
        harmonic_solution = df;

        % Create a logical mask for the missing values
        missing_points = ~known_points;

        % Compute weights for known points
        weights = abs(max(curr_param(:, :, 5) + curr_param(:, :, 6), ...
                          curr_param(:, :, 5) - curr_param(:, :, 6)));
        weights(~known_points) = 0; % Set weights to zero for missing points

        % Iteratively solve the Laplace equation for the missing points
        tolerance = 1e-6;
        max_iterations = 1000;
        iteration = 0;
        diff = Inf;

        while diff > tolerance && iteration < max_iterations
            old_solution = harmonic_solution;

            % Update only the missing points using the weighted Laplace stencil
            for i = 2:size(harmonic_solution, 1) - 1
                for j = 2:size(harmonic_solution, 2) - 1
                    if missing_points(i, j)
                        neighbor_avg = (harmonic_solution(i+1, j) + ...
                                        harmonic_solution(i-1, j) + ...
                                        harmonic_solution(i, j+1) + ...
                                        harmonic_solution(i, j-1)) / 4;
                        
                        % Incorporate the weighted contribution of known values
                        harmonic_solution(i, j) = ...
                            (weights(i, j) * df(i, j) + neighbor_avg) / ...
                            (weights(i, j) + 1);
                    end
                end
            end

            % Compute the maximum change between iterations
            diff = max(abs(harmonic_solution(:) - old_solution(:)));
            iteration = iteration + 1;
        end

        % Store the interpolated slice
        pdf(:, :, index) = harmonic_solution;
    end

    % Ensure the mask is respected and original data remains in masked regions
    for x = 1:size(mask, 1)
        for y = 1:size(mask, 2)
            if ~mask(y, x)
                pdf(y, x, :) = curr_param(y, x, :);
            end
        end
    end

    % Update the parameters field with the interpolated data
    function_outputs.parameters = pdf;
end
