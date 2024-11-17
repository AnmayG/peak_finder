function function_outputs = spline(function_inputs)
    % Initialize function_outputs as the first input struct
    function_outputs = function_inputs(1);
    parameters_dataframe = function_inputs(1).parameters;
    in = function_inputs(1);
    mask = in.new_settings(1).roi_mask;
    curr_param = function_inputs(1).parameters;

    % Iterate over each slice in the third dimension
    for index = 1:size(parameters_dataframe, 3)
        % Extract the current slice
        df = parameters_dataframe(:, :, index);

        % Find known points
        [x, y] = meshgrid(1:size(df, 2), 1:size(df, 1));
        known_points = df ~= 0;  % Known points (non-zero values)
        
        % If no known points, skip this slice
        if ~any(known_points, 'all')
            continue;
        end

        % Flatten known points to vectors
        x_known = x(known_points);
        y_known = y(known_points);
        values_known = df(known_points);

        % Perform cubic spline interpolation
        interpolated_df = griddata(x_known(:), y_known(:), values_known(:), x, y, 'cubic');

        % Replace the original slice with the interpolated data
        pdf(:, :, index) = interpolated_df;
    end

    %disp(size(mask))
    for x=1:size(mask,1)
        for y=1:size(mask,2)
           if(~mask(y, x))
                pdf(y, x, :) = curr_param(y, x, :);
                
           end
        end
    end

    % Update the parameters field with the interpolated data
    function_outputs.parameters = pdf;
end
