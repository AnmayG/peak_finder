% Assuming zz, xx, and yy are already defined and have the same dimensions.
zz=data;
zzz=zz-1;
[xxx,yyy]=MITBxByFromBz(zzz,236);
% Create the figure
figure;

magnitude = sqrt(xxx.^2 + yyy.^2 + zzz.^2);

% Define the threshold for vector length
threshold = 0.1; % Adjust this value as needed

% Filter vectors based on their magnitude
valid = magnitude <= threshold; % Logical array of valid vectors
xxx(~valid) = 0; % Set invalid vectors to zero
yyy(~valid) = 0;
%zzz(~valid) = 0;

% Plot the heatmap
imagesc(zzz); % Display heatmap of zz
clim([-0.2 0.2]);
colormap winter; % Choose a colormap for the heatmap
colorbar;    % Add a colorbar for reference
axis xy;     % Ensure correct orientation
axis off;    % Turn off the axis
set(gcf, 'Color', 'w'); % Set the background color to white

% Overlay the vector field
hold on;

% Define the grid for quiver
[X, Y] = meshgrid(1:size(zz, 2), 1:size(zz, 1));

% Plot every 'step' vector for clarity
step = 10; % Adjust this for sparser or denser arrows
quiver(X(1:step:end, 1:step:end), Y(1:step:end, 1:step:end), ...
       xxx(1:step:end, 1:step:end), yyy(1:step:end, 1:step:end), ...
       'w', 'AutoScale', 'on', 'AutoScaleFactor', 1); % 'w' sets arrow color to white

hold off;

% Adjust figure appearance
axis square; % Make the axis square for a clean look
title('Heatmap with Vector Field Overlay');
%% 
% Assuming xx, yy, zz are 3D vector components defined over a grid.

% Define the grid (X, Y, Z)
[X, Y, Z] = meshgrid(1:size(xx, 2), 1:size(xx, 1), 1:size(xx, 3));
% zz=data;
% zzz=zz-1;
% Calculate the magnitude of the vectors
magnitude = sqrt(xxx.^2 + yyy.^2 + zzz.^2);

% Define the threshold for vector length
threshold = 0.1; % Adjust this value as needed

% Filter vectors based on their magnitude
valid = magnitude <= threshold; % Logical array of valid vectors
xxx(~valid) = 0; % Set invalid vectors to zero
yyy(~valid) = 0;
zzz(~valid) = 0;

% Create the figure
figure;
hold on;

% Plot the 3D vector field (only valid vectors are shown)
step = 3; % Adjust step size for sparser vectors
quiver3(X(1:step:end, 1:step:end, 1:step:end), ...
        Y(1:step:end, 1:step:end, 1:step:end), ...
        Z(1:step:end, 1:step:end, 1:step:end), ...
        xx(1:step:end, 1:step:end, 1:step:end), ...
        yy(1:step:end, 1:step:end, 1:step:end), ...
        zzz(1:step:end, 1:step:end, 1:step:end), ...
        'AutoScale', 'on', 'AutoScaleFactor', 10, 'Color', 'b'); % Adjust color and scale

% Customize the axes and view
axis tight;
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Vector Field with Threshold Filtering');

% Set background color and grid
set(gcf, 'Color', 'w'); % Set the figure background to white
grid on;
view(3); % Set to 3D view
hold off;
