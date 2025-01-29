% Get the dimensions of f2 parameters
[~, shape_x, shape_y] = size(f2.parameters);

% Extract p2 for mode calculation
original_p2 = round(f2.parameters(2, :, :), 3);

% Flatten p2 to calculate the mode
original_p2_flat = reshape(original_p2, [], 1);
mode_value = mode(original_p2_flat); % Calculate the mode of the rounded values

% Initialize the output array for plot values
plot_values = zeros(shape_x, shape_y);

% Loop through all x, y coordinates
for x = 1:shape_x
    for y = 1:shape_y
        % Extract parameters for clarity
        p2 = f2.parameters(2, x, y);
        p3 = f2.parameters(3, x, y);

        % Compute the three possible values
        p2_no_change = p2;
        p2_add = p2 + p3;
        p2_subtract = p2 - p3;

        % Round to 3 decimal places for comparison
        p2_no_change_rounded = round(p2_no_change, 3);
        p2_add_rounded = round(p2_add, 3);
        p2_subtract_rounded = round(p2_subtract, 3);

        % Compute the absolute differences with the mode value
        diff_no_change = abs(p2_no_change_rounded - mode_value);
        diff_add = abs(p2_add_rounded - mode_value);
        diff_subtract = abs(p2_subtract_rounded - mode_value);

        % Find the action that minimizes the difference
        [~, idx] = min([diff_no_change, diff_add, diff_subtract]);

        % Apply the corresponding adjustment
        if idx == 1
            plot_values(x, y) = p2_no_change;
        elseif idx == 2
            plot_values(x, y) = p2_add;
        else
            plot_values(x, y) = p2_subtract;
        end
    end
end

% Plot the results
imagesc(plot_values'); 
minValue = 0.571; 
maxValue = 0.581; 
caxis([minValue maxValue]); 
colorbar;

%colormap('viridis'); % Change colormap as needed
% title('Adjusted Plot to Minimize Difference to Mode');
% xlabel('y');
% ylabel('x');
% set(gca, 'YDir', 'normal'); % Ensure the y-axis is not flipped
