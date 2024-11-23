function function_outputs = refit_red(function_inputs)
    %   function_inputs = struct array with fields:
    %       signal: [n × n × freq double]
    %   SweepParam: freq
    %   parameters: [n × n × 17 double]
    %    peak_info: [1×1 struct] with fields locs, vals, widths, proms
    %   reflection: [n × n double]
    %     settings: [1×1 struct] with each parameter from inputs
    % new_settings: [1×1 struct] with parameters at time of run
    function_outputs = function_inputs;
    settings = function_inputs.new_settings;
    % Modify function_outputs fields below
    num_inputs = length(function_inputs);
    for i=1:num_inputs
        in = function_inputs(i);
        mask = in.new_settings(i).roi_mask;
        prev_params = in.parameters;
        for x=1:size(mask,1)
            for y=1:size(mask,2)
                if mask(y,x) && prev_params(y, x, 1) ~= 0 % D map is nonzero
                    mask(y,x) = 0;
                end
            end
        end
        [pdf, ~, peak_info] = full_peak_wrapper([], in.signal, in.SweepParam, ...
            settings, 0, in.reflection, 16);
        curr_param = function_inputs(i).parameters;
        for x=1:size(mask,1)
            for y=1:size(mask,2)
               if(~mask(y, x))
                    pdf(y, x, :) = curr_param(y, x, :);    
               end
            end
        end
        function_outputs(i).parameters = pdf;
        function_outputs(i).peak_info = peak_info;
    end
end