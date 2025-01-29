function function_outputs = vector_B(function_inputs)
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
        em=in.parameters(:,:,2);
        bzm=em/2.8e6;
        % Get the bounding box of the largest rectangular region within mask
        [y_indices, x_indices] = find(mask);
        if isempty(y_indices) || isempty(x_indices)
            error('Mask is empty, cannot extract a rectangular subregion.');
        end

        % Determine bounding box
        y_min = min(y_indices);
        y_max = max(y_indices);
        x_min = min(x_indices);
        x_max = max(x_indices);

        % Extract rectangular subregion from relevant fields
        bz = bzm(y_min:y_max, x_min:x_max, :);
        assignin('base', 'bz', bz);
        fs = min(size(bz, 1), size(bz, 2));
        [bx,by]=MITBxByFromBz(bz,fs);


        figure;

        magnitude = sqrt(bx.^2 + by.^2 + bz.^2);
        
        % Define the threshold for vector length
        threshold = 0.1; % Adjust this value as needed
        
        % Filter vectors based on their magnitude
        valid = magnitude <= threshold; % Logical array of valid vectors
        % bx(~valid) = 0; % Set invalid vectors to zero
        % by(~valid) = 0;
        %zzz(~valid) = 0;
        assignin('base', 'bx', bx);
        assignin('base', 'by', by);
        % Plot the heatmap
        imagesc(bz); % Display heatmap of zz
        %clim([-0.2 0.2]);
        colormap winter; % Choose a colormap for the heatmap
        colorbar;    % Add a colorbar for reference
        axis xy;     % Ensure correct orientation
        axis off;    % Turn off the axis
        set(gcf, 'Color', 'w'); % Set the background color to white
        
        % Overlay the vector field
        hold on;
        
        % Define the grid for quiver
        [X, Y] = meshgrid(1:size(bz, 2), 1:size(bz, 1));
        
        % Plot every 'step' vector for clarity
        step = 10; % Adjust this for sparser or denser arrows
        quiver(X(1:step:end, 1:step:end), Y(1:step:end, 1:step:end), ...
               bx(1:step:end, 1:step:end), by(1:step:end, 1:step:end), ...
               'w', 'AutoScale', 'on', 'AutoScaleFactor', 1); % 'w' sets arrow color to white
        
        hold off;
        
        % Adjust figure appearance
        axis square; % Make the axis square for a clean look
        title('Heatmap with Vector Field Overlay');


    end
end