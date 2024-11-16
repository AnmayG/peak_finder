function function_outputs = roi_to_red(function_inputs)
    %   function_inputs = struct array with fields:
    %       signal: [n × n × freq double]
    %   SweepParam: freq
    %   parameters: [n × n × 17 double]
    %    peak_info: [1×1 struct] with fields locs, vals, widths, proms
    %   reflection: [n × n double]
    %     settings: [1×1 struct] with each parameter from inputs
    % new_settings: [1×1 struct] with parameters at time of run
    function_outputs = function_inputs;
    % Modify function_outputs fields below
    num_inputs = length(function_inputs);
    mask = function_inputs.new_settings.roi_mask;
    for i=1:num_inputs
        pdf = function_outputs(i).parameters;
        red_out = zeros(size(pdf));
        for x=1:size(mask,1)
            for y=1:size(mask,2)
                if(mask(y,x))
                    pdf(y,x,:) = red_out(y,x,:);
                end
            end
        end
        function_outputs(i).parameters = pdf;
    end
end
