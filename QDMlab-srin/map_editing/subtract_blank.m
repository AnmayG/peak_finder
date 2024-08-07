function subtractedData = subtract_blank(kwargs)
%[subtractedData] = subtract_blank('nFiles', 'blankFile', 'alignmentType', 'checkPlot', 'save')
% Subtracts a blank map from the Data
%
% Parameters
% ----------
%   nFiles: ['none']
%       A path (char) or cell of paths that contain the data and are to be
%       corrected.
%   blankFile: path ['none']
%       The path of the blank B111 map.
%   checkPlot: bool [false]
%       if true: creates a plot to check if the subtraction worked as
%       expected
%   save: bool [true]
%       if true: a new file with the name: `B111BlankSub.mat` will be
%       created. Otherwise the results will only be returned.
%   alignmentType: str ['laser']
%       laser: laser images are used for alignment -> subtraction of blank maps
%       led: led images are used for alignment -> subtraction of measurements
%
% Note
% ----
%   If no :code:`nFiles` or :code:`blankFile` is passed to the function,
%   you will be prompted to select them.

arguments
    kwargs.nFiles = 'none'
    kwargs.blankFile = 'none'
    kwargs.alignmentType = 'laser'
    kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)} = false
    kwargs.save {mustBeBoolean(kwargs.save)} = true
end

subtractedData = containers.Map();

laserFileName = 'laser.jpg';

%% manual
% if nFoilders and blankData uses default values i.e. false
nFiles = automatic_input_ui__(kwargs.nFiles, 'title', 'Select measurement file', 'type', 'file');
blankFile = automatic_input_ui__(kwargs.blankFile, 'title', 'Select blank file', 'single', true, 'type', 'file');

%% automatic subtraction for all folders
% checks if none of the default arguments is used
nFiles = correct_cell_shape(nFiles);

msg = sprintf('loading blank file: << %s >>', blankFile);
logMsg('info',msg,1,0);

blankData = load(blankFile);

% [nTransForms, nRefFrames] = get_tform_multi(blankFolder, nFolders, ...
%                             'reverse', true, ...
%                             'laser',true, 'checkPlot', checkPlot);
[blankFolder ]= fileparts(blankFile);

[bool, dataName,ledName] = is_B111(blankData);

switch kwargs.alignmentType
case 'laser'
    if isfield(blankData, 'laser')
        movingData = blankData.laser;
    else
        movingData = imread(fullfile(blankFolder, laserFileName));
    end
case 'led'
    movingData = blankData.(ledName);
end

for i = 1 : size(nFiles, 2)
    iFile = nFiles{i};
    [iFolder fileName ext] = fileparts(iFile);

    msg = sprintf('loading laser & magnetic data: << %s >>', iFolder);
    logMsg('info',msg,1,0);
    
    fileData = load(iFile);

    switch kwargs.alignmentType
    case 'laser'
        if isfield(fileData, 'laser')
            fixedData = fileData.laser;
        else
            fixedData = imread(fullfile(iFolder, laserFileName));
        end
    case 'led'
        fixedData = fileData.(ledName);
    end

    % fixedData = imread(fullfile(iFolder, laserFileName));
    % copy data for saving
    newFileData = fileData;
    
    [transForm, refFrame] = get_image_tform(fixedData, movingData,...
        'checkPlot', kwargs.checkPlot, 'title', 'laser alignment');

    % transform the blank
    B111ferroTransformed = tform_data(blankData.B111ferro, transForm, refFrame);
    B111paraTransformed = tform_data(blankData.B111para, transForm, refFrame);

    % crop the FOV and subtract blank
    % [x, y, w, h] = get_mask_extent(B111ferroTransformed);
    % B111ferro = fileData.B111ferro;
    
    B111ferroTransformed(B111ferroTransformed==0) = nan;
    fileB111ferro = fileData.B111ferro;
    fileB111para = fileData.B111para;
    
    newFileData.B111ferro = fileB111ferro- B111ferroTransformed;
    newFileData.B111para = fileB111para - B111paraTransformed;
%     B111para = fileData.B111para - B111paraTransformed;
    newFileData.blank = blankFile;
    newFileData.blankB111ferro = B111ferroTransformed;
    newFileData.blankB111Para = B111paraTransformed;
    subtractedData(iFolder) = newFileData;
    
    if kwargs.save
        saveFilePath = fullfile(iFolder, sprintf('%s-blank.mat', fileName));
        msg = sprintf('%s', saveFilePath);
        logMsg('SAVE',msg,1,0);
        save(saveFilePath, '-struct', 'newFileData');
    end

    if kwargs.checkPlot
        med = abs(median(fileB111ferro, 'all', 'omitnan'));
        st = std(fileB111ferro, [], 'all', 'omitnan');
        
        fig = figure('Units', 'normalized', ...
                 'Position',[0.1 0.1 0.8 0.8], 'Name', 'Blank Subtraction');
        
        sp1 = subplot(2,2,1);
        QDM_figure(fileB111ferro, 'ax', sp1, 'title','Original');
        sp2 = subplot(2,2,2);
        QDM_figure(blankData.B111ferro, 'ax', sp2, 'title','Blank');        
        sp3 = subplot(2,2,3);
        QDM_figure(fileData.B111ferro - blankData.B111ferro, 'ax', sp3, 'title','Unaligned subtraction');     
        sp4 = subplot(2,2,4);
        QDM_figure(newFileData.B111ferro, 'ax', sp4, 'title','Aligned subtraction', 'preThreshold','none');

        ax = [sp1, sp2, sp3, sp4];
        linkaxes(ax)
        
        for a = ax
           colorbar(a)
%            set(ax, 'CLim', [-1, 1]*(med + 4 * st));
        end
    end
end
