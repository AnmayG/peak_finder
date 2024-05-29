load('..\LaH10\LaH10_-23dB_500mW_146K_50G.mat');

split = 4.21e9;
tolerance = 1.5e8;
xsize = 80;
ysize = 80;
num_params = 13;
parameters_dataframe = zeros(ysize, xsize, num_params);
for x = 1:xsize
    for y = 1:ysize
        disp(x + " " + y)
        seed = peak_find_function(split, tolerance, x, y, gWide);
        parameters_dataframe(y, x, :) = seed;
    end
end

fprintf("Completed, generating figure\n")

parameter_labels = ["D_{111} (Hz)", "E_{111} (Hz)", "Width (σ)", "111 Width Deviation (σ)", "111 Contrast (%)", "111 Contrast Deviation", ... 
"D_{non111} (Hz)", "E_{non111}", "Non-111 Width (σ)", "Non-111 Width Deviation", "Non-111 Contrast (%)", "Non-111 Contrast Deviation"];

% Create a figure
figure;

% Loop through each slice and display it
for i = 1:12
    subplot(3, 4, i); % Create a subplot grid of 3 rows and 4 columns
    data = parameters_dataframe(:, :, i);
    [zeroRows, zeroCols] = find(data == 0);
    data(data == 0) = NaN;
    imagesc(data); % Display the slice
    colorbar; % Add a colorbar to each subplot
    hold on;
    plot(zeroCols, zeroRows, 'rs', 'MarkerSize', 1, 'MarkerFaceColor', 'r');
    % Hold off the plot
    hold off;
    title([parameter_labels(i)]); % Add a title to each subplot
    axis equal; % Set axis scales to be equal,
    axis tight; % Remove excess white space around the image
end

e111_50g_146k = squeeze(parameters_dataframe(:, :, 2)) .^ 2;

load('..\LaH10\LaH10_-23dB_500mW_146K_0G.mat');

split = 4.21e9;
tolerance = 1.5e8;
xsize = 80;
ysize = 80;
num_params = 13;
parameters_dataframe = zeros(ysize, xsize, num_params);
for x = 1:xsize
    for y = 1:ysize
        disp(x + " " + y)
        seed = peak_find_function(split, tolerance, x, y, gWide);
        parameters_dataframe(y, x, :) = seed;
    end
end

fprintf("Completed, generating figure\n")

parameter_labels = ["D_{111} (Hz)", "E_{111} (Hz)", "Width (σ)", "111 Width Deviation (σ)", "111 Contrast (%)", "111 Contrast Deviation", ... 
"D_{non111} (Hz)", "E_{non111}", "Non-111 Width (σ)", "Non-111 Width Deviation", "Non-111 Contrast (%)", "Non-111 Contrast Deviation"];

% Create a figure
figure;

% Loop through each slice and display it
for i = 1:12
    subplot(3, 4, i); % Create a subplot grid of 3 rows and 4 columns
    data = parameters_dataframe(:, :, i);
    [zeroRows, zeroCols] = find(data == 0);
    data(data == 0) = NaN;
    imagesc(data); % Display the slice
    colorbar; % Add a colorbar to each subplot
    hold on;
    plot(zeroCols, zeroRows, 'rs', 'MarkerSize', 1, 'MarkerFaceColor', 'r');
    % Hold off the plot
    hold off;
    title([parameter_labels(i)]); % Add a title to each subplot
    axis equal; % Set axis scales to be equal,
    axis tight; % Remove excess white space around the image
end

e111_0g_146k = squeeze(parameters_dataframe(:, :, 2)) .^ 2;

data = sqrt(e111_50g_146k - e111_0g_146k) ./ 1.4e8;
data(data ~= real(data)) = 6;
figure;
hold on;
title("Quadrature Difference between 50G and 0G at 146K")
imagesc(data)
axis equal;
axis tight;
colorbar;
clim([0, 1.3])
hold off;