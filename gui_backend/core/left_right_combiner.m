function gFile = left_right_combiner(f1, f2, filePath)
    % COMBINE_GWIDE_FILES Combines two gWide structures with linear interpolation.
    % 
    % gWide = combine_gWide_files(leftFile, rightFile)
    % 
    % Inputs:
    %   leftFile  - Filename (including path) of the left .mat file
    %   rightFile - Filename (including path) of the right .mat file
    %
    % Output:
    %   gWide - Combined gWide structure with interpolated signal
    
    % Load the left and right .mat files

    % Determine which file contains "left" and which contains "right"
    if contains(f1, 'left') && contains(f2, 'right')
        leftFile = fullfile(filePath, f1);
        rightFile = fullfile(filePath, f2);
    elseif contains(f1, 'right') && contains(f2, 'left')
        leftFile = fullfile(filePath, f2);
        rightFile = fullfile(filePath, f1);
    else
        error('Selected files must contain "left" and "right" in their names.');
    end
    leftData = load(leftFile);
    rightData = load(rightFile);
    
    % Access the gWide structure from both files
    gWide_left = leftData.gWide;
    gWide_right = rightData.gWide;
    
    % Normalize the signal field of both gWide_left and gWide_right
    % Normalization is based on the mean of the first 5 elements of each pixel's z-dimension
    
    % Normalize gWide_left.signal
    mean_first_5_left = mean(gWide_left.signal(:,:,1:5), 3);
    gWide_left.signal = gWide_left.signal ./ mean_first_5_left;
    
    % Normalize gWide_right.signal
    mean_first_5_right = mean(gWide_right.signal(:,:,1:5), 3);
    gWide_right.signal = gWide_right.signal ./ mean_first_5_right;
    
    % Calculate the z-spacing from gWide_left's SweepParam
    left_z_spacing = mean(diff(gWide_left.SweepParam));
    
    % Determine the required number of interpolation steps
    z_gap = gWide_right.SweepParam(1) - gWide_left.SweepParam(end);
    num_interp_steps = round(z_gap / left_z_spacing);
    
    % Initialize the interpolated signal array
    interp_signal = zeros(size(gWide_left.signal, 1), size(gWide_left.signal, 2), num_interp_steps);
    
    % Parallel interpolation calculation
    parfor z = 1:num_interp_steps
        % Compute weights for linear interpolation
        weight_left = (num_interp_steps - z + 1) / (num_interp_steps + 1);
        weight_right = z / (num_interp_steps + 1);
        % Perform weighted average for each z-slice in the interpolation range
        interp_signal(:,:,z) = gWide_left.signal(:,:,end) * weight_left + gWide_right.signal(:,:,1) * weight_right;
    end
    
    % Append the SweepParam, signal fields, and interpolated gap
    interp_sweep_param = linspace(gWide_left.SweepParam(end), gWide_right.SweepParam(1), num_interp_steps);
    gWide_left.SweepParam = [gWide_left.SweepParam, interp_sweep_param, gWide_right.SweepParam];
    gWide_left.signal = cat(3, gWide_left.signal, interp_signal, gWide_right.signal);
    gWide = gWide_left;
    % Return the combined gWide structure
    % Create the new filename by replacing 'left' with 'cbnd' in the left filename
    [~, leftFileName, ~] = fileparts(leftFile);  % Get the left file's name without extension
    newFileName = strrep(leftFileName, 'left', 'cbnd');  % Replace 'left' with 'cbnd'
    newFileFullPath = fullfile(filePath, [newFileName, '.mat']);  % Add .mat extension
    % disp(newFileFullPath)
    % assignin('base', 'huh', newFileFullPath{1});
    % Save the modified gWide structure into the new file
    disp(newFileFullPath)
    save(newFileFullPath, 'gWide');
    
    disp(['File saved as: ', newFileFullPath]);

    gFile=newFileFullPath;
end
