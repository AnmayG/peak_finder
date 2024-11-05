function function_outputs = gaussian_blur(function_inputs)
    % function_inputs is a struct array with fields:
    %     signal: [n × n × freq double]
    % SweepParam: freq
    % parameters: [n × n × 17 double]
    %  peak_info: [1×1 struct] with fields locs, vals, widths, proms
    % reflection: [n × n double]
    %   settings: [1×1 struct] with each parameter from inputs
    % function_outputs should also be in this format, but only 1 struct

    function_outputs = function_inputs(1);
    parameters_dataframe = function_inputs(1).parameters;

    % Iterate over each "slice" in the third dimension
    for index = 1:size(parameters_dataframe, 3)
        % Extract the current slice
        df = parameters_dataframe(:, :, index);

        % Find known and unknown points
        [x, y] = meshgrid(1:size(df, 2), 1:size(df, 1));
        known_points = df ~= 0; % Known points (non-zero values)
        x_known = x(known_points);
        y_known = y(known_points);
        values_known = df(known_points);

        % Perform cubic spline interpolation
        interpolated_df = griddata(x_known, y_known, values_known, x, y, 'cubic');

        % Replace the original slice with the interpolated data
        parameters_dataframe(:, :, index) = interpolated_df;
    end

    % Update the parameters field with the interpolated data
    function_outputs.parameters = parameters_dataframe;
end
