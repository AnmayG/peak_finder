load('..\327_20gpa_lowT_vsB\327_20gpa_20K_4A.mat');

split = 3.1e9;
tolerance = 1.0e8;
xsize = 180;
ysize = 180;
num_params = 16;
parameters_dataframe = zeros(ysize, xsize, num_params);
data = gWide.signal;
freq = gWide.SweepParam;
for x = 1:xsize
    for y = 1:ysize
        disp(x + " " + y)
        seed = peak_find_function(split, tolerance, x, y, data, freq);
        parameters_dataframe(y, x, :) = seed;
    end
end

laplacian_smoothing(parameters_dataframe, true)

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
