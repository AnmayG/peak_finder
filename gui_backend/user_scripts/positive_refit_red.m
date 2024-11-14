function function_outputs = positive_refit_red(function_inputs)
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
    % assignin('base', 'setz', settings);
    % Modify function_outputs fields below
    num_inputs = length(function_inputs);
    for i=1:num_inputs
        in = function_inputs(i);
        mask = in.settings(i).roi_mask;


        prev_params = in.parameters;
        for x=1:size(mask,1)
            for y=1:size(mask,2)
                if mask(y,x) && prev_params(y, x, 1) ~= 0 % D map is nonzero
                    mask(y,x) = 0;
                end
            end
        end
        sig=in.signal;
        %no longer flip here
        % mn = median(sig, 3);
        % 
        % % Replicate mn along the third dimension to match the dimensions of signal
        % mn_expanded = repmat(mn, 1, 1, size(sig, 3));
        % 
        % % Replace signal(x, y, z) with 2 * mn - signal(x, y, z)
        % flip_signal = 2 * mn_expanded - sig;
        % 
        flip_signal=sig;


        [pdf, ~, peak_info] = full_peak_wrapper([], flip_signal, in.SweepParam, ...
            settings, 0, in.reflection, 16);

        % assignin('base', 'sig', flip_signal);
        % %disp('hi');
        % assignin('base', 'sp', in.SweepParam);
        curr_param = function_inputs(i).parameters;
        for x=1:size(mask,1)
            for y=1:size(mask,2)
               if(~mask(y, x))
                    pdf(y, x, :) = curr_param(y, x, :);
               else
                   pdf(y, x, 2) = 2.863*pdf(y, x, 2); %2.863?
               end
            end
        end
        function_outputs(i).parameters = pdf;
        % assignin('base', 'pom', function_outputs(i).parameters);


        function_outputs(i).peak_info = peak_info;
    end
end