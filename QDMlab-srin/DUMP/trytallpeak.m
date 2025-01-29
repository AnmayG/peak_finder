% Get the dimensions of f2 parameters

[~, shape_x, shape_y] = size(f2.parameters);

% Initialize the output array for plot values
plot_values = zeros(shape_x, shape_y);

% Loop through all x, y coordinates
for x = 1:shape_x
    for y = 1:shape_y
        % Extract parameters for clarity
        p2 = f2.parameters(2, x, y);
        p3 = f2.parameters(3, x, y);
        p4 = f2.parameters(4, x, y);
        p5 = f2.parameters(5, x, y);
        p6 = f2.parameters(6, x, y);
        p7 = f2.parameters(7, x, y);

        % Compute the conditions
        condition1 = p4 - p5;
        condition2 = p6 - p7;

        % Define X (p4 and p5 within a factor of 2)
        r=0.33;
        X = (abs(p4) <= 1/r * abs(p5) && abs(p4) >= r * abs(p5));

        % Define Y (p6 and p7 within a factor of 2)
        Y = (abs(p6) <= 1/r * abs(p7) && abs(p6) >= r * abs(p7));

        % Apply the logic based on the conditions
        if sign(condition1) == sign(condition2) || (X && Y) % formerly &&
            plot_values(x, y) = p2;
        elseif condition1 > 0 && condition2 < 0
            plot_values(x, y) = p2 - p3;
        else
            plot_values(x, y) = p2 + p3;
        end
    end
end

% Plot the results
imagesc(plot_values'); minValue = 0.571; maxValue = 0.581; caxis([minValue maxValue]); colorbar;

%colormap('viridis'); % Change colormap as needed
% title('Plot based on conditions');
% xlabel('y');
% ylabel('x');
% set(gca, 'YDir', 'normal'); % Ensure the y-axis is not flipped
