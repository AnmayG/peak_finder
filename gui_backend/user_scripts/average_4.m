function function_outputs = average_4(function_inputs)
    %   function_inputs = struct array with fields:
    %       signal: [n × n × freq double]
    %   SweepParam: freq
    %   parameters: [n × n × 17 double]
    %    peak_info: [1×1 struct] with fields locs, vals, widths, proms
    %   reflection: [n × n double]
    %     settings: [1×1 struct] with each parameter from inputs
    % function_outputs should also be in this format, but only 1 struct
    k_rad = 4;
    function_outputs = function_inputs(1);
    K = (1/k_rad^2)*ones(k_rad);
    parameters_dataframe = function_inputs(1).parameters;
    for index = 1:size(parameters_dataframe, 3)
        parameters_dataframe(:, :, index) = conv2(parameters_dataframe(:, :, index), K, 'same');
    end
    function_outputs.parameters = parameters_dataframe;
end