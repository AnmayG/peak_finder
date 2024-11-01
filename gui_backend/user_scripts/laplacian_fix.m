function function_outputs = laplacian_fix(function_inputs)
    %   function_inputs = struct array with fields:
    %       signal: [n × n × freq double]
    %   SweepParam: freq
    %   parameters: [n × n × 17 double]
    %    peak_info: [1×1 struct] with fields locs, vals, widths, proms
    %   reflection: [n × n double]
    %     settings: [1×1 struct] with each parameter from inputs
    % function_outputs should also be in this format, but only 1 struct
    function_outputs = function_inputs(1);
    [pdf, ~] = laplacian_wrapper(function_outputs.settings, ...
        function_outputs.parameters, function_outputs.peak_info, 0);
    % Plot difference to check if it actually does something
    % figure;
    % imagesc(squeeze(function_inputs(1).parameters(:,:,1) == pdf(:,:,1)));
    function_outputs.parameters = pdf;
end