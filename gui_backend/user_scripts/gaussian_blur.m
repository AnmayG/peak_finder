function function_outputs = gaussian_blur(function_inputs)
    %   function_inputs = struct array with fields:
    %       signal: [n × n × freq double]
    %   SweepParam: freq
    %   parameters: [n × n × 17 double]
    %    peak_info: [1×1 struct] with fields locs, vals, widths, proms
    %   reflection: [n × n double]
    %     settings: [1×1 struct] with each parameter from inputs
    % function_outputs should also be in this format, but only 1 struct
    function_outputs = function_inputs(1);
    parameters_dataframe = function_inputs(1).parameters;
    for index = 1:size(parameters_dataframe, 3)
        df = parameters_dataframe(:, :, index);
        df(isnan(df)) = 0;
        parameters_dataframe(:, :, index) = imgaussfilt(df, 3);
    end
    function_outputs.parameters = parameters_dataframe;
end