function function_outputs = queue_fits(function_inputs)
    %   function_inputs = struct array with fields:
    %       signal: [n × n × freq double]
    %   SweepParam: freq
    %   parameters: [n × n × 17 double]
    %    peak_info: [1×1 struct] with fields locs, vals, widths, proms
    %   reflection: [n × n double]
    %     settings: [1×1 struct] with each parameter from inputs
    % new_settings: [1×1 struct] with parameters at time of run
    function_outputs = function_inputs;
    params = function_inputs.new_settings;
    % Modify function_outputs fields below
    num_inputs = length(function_inputs);
    for i=1:num_inputs
        in = function_inputs(i);
        [pdf, ~, peak_info] = full_peak_wrapper([], in.signal, in.SweepParam, ...
            params, 0, in.reflection, 16);
        function_outputs(i).parameters = pdf;
        function_outputs(i).peak_info = peak_info;
    end
end