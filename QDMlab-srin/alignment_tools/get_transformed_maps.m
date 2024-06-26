function [transformedData, nFiles] = get_transformed_maps(nFolders, kwargs, filter)
%[transformedData, nFiles] = get_transformed_maps(nFolders; 'fileName', 'transFormFile', 'fixedIdx', 'upCont', 'quadBgSub', 'reverse', 'checkPlot', 'filterProps')
% These codes (1) register the maps and (2) analizes a user selected magnetic
% pattern for changes from one map to the next.(folders, varargin)
%
% Parameters
% ----------
%     nFolders:list
%         list of absolute path to each data folder. First entry in the list
%         is used as reference image.
%         Note: calculations in the order the folders are given
%     kwargs.fileName: str ['Bz_uc0']
%         name of the .mat file to be used.
%     kwargs.transFormFile: str ['none']
%         absolute path of a file where previously calculated transformations
%         are in. If file does not exist the file will be created so that
%         you dont have to do this every time.
%     removeHotPixel: double [0]
%         Uses 'filter_hot_pixel' function. All filtering is done pre
%         tranfornation.
%
%         if false: hotpixels will not be removed
%         else:     removeHotPixel determines how many standard deviations
%                   have to be exceeded for the pixel to be filtered
%     includeHotPixel: bool [0]
%           uses 'filter_hot_pixel' function
%         if true: hotpixels is included in the calculation of the mean
%         if false: hotpixels is NOT included in the calculation of the mean
%     winSize: int [4]
%         Window size around the hot pixel for mean calculation (e.g. 4 ->
%         4x4 window)
%     chi: bool [0]
%         if true: chi^2 values are used to filter (sum of target chi^2)
%         if false: data is filtered by data values
%     kwargs.fixedIdx: int [1]
%         index of the reference LED image. This will fixed while the other
%         image is transformed.
%     reverse: bool [0]
%         NOT IMPLEMENTED YET
%         if true:  refernce - tform -> target
%         if false: target   - tform -> reference
%
% See also, filter_hot_pixels

arguments
    nFolders;
    kwargs.fileName = 'Bz_uc0';
    kwargs.transFormFile = 'none';
    kwargs.fixedIdx = 1;
    kwargs.upCont = false;
    kwargs.quadBgSub (1,1) {mustBeBoolean(kwargs.quadBgSub)} = false;

    % other
    kwargs.reverse = false;
    kwargs.checkPlot = false;
    % filter related
    filter.filterProps struct = struct();
end

nFolders = correct_cell_shape(nFolders);

%% tranformation / filtering
transformedData = containers.Map;
nFiles = {};

% define QDM parameters
pixelsize = 4.68e-6;

% check if filename has '.mat'
kwargs.fileName = check_suffix(kwargs.fileName);

% generate reference file name
refFile = fullfile(nFolders{kwargs.fixedIdx}, filesep, kwargs.fileName);
refFileData = load(refFile);

% get transformations and rframes
if kwargs.checkPlot
    [nTransForms, nRefFrames] = align_images(nFolders, 'transformFile', 0, ...
        'fileName', kwargs.fileName, 'fixedIdx', kwargs.fixedIdx);
else
    [nTransForms, nRefFrames] = get_tform_multi(refFile, nFolders, ...
        'transFormFile', kwargs.transFormFile, 'reverse', kwargs.reverse);
end

[bool, dataName,ledName] = is_B111(refFileData);
refData = refFileData.(dataName);
refLed = refFileData.(ledName);

% cycle through all folders
for i = 1:size(nFolders, 2)
    % create filename
    iFolder = nFolders{i};
    iFile = fullfile(iFolder, filesep, kwargs.fileName);

    msg = sprintf('loading << %s >> target file for transformation', iFile');
    logMsg('debug', msg, 1, 0);
    %     fprintf('<> loading << %s >> target file for transformation\n', iFile(end-size(iFile,2)/2:end))

    nFiles{end+1} = iFile;

    target = load(iFile);
    
    [~, dataName, ledName] = is_B111(target);
    % read data and threshold to 5
    targetData = target.(dataName);
    targetLed = target.(ledName);
    
    if kwargs.quadBgSub
        targetData = QuadBGsub(targetData);
    end
    % pre filtering
    targetData = filter_hot_pixels(targetData);

    %% upward cont.
    if iscell(kwargs.upCont)
        h = kwargs.upCont{i};
        if h ~= 0
            msg = sprintf('calculating upward continuation (%i) micron', h);
            logMsg('info', msg, 1, 0);
            %             disp(['<>   calculating upward continuation (' num2str(h) ') micron'])
            targetData = UpCont(targetData, h*1e-6, 1/pixelsize);
        end
    end

    %% filtering
    if ~all(structfun(@isempty, filter.filterProps))
        if isfield(filter.filterProps, 'chi') && isequal(filter.filterProps.chi, true)
            filter.filterProps.chi = target.chi2Pos1 + target.chi2Pos2 + target.chi2Neg1 + target.chi2Neg2;
        end
        filterProps = namedargs2cell(filter.filterProps);
        msg = sprintf('filtering: ...%s... .mat', iFile(end -40:end - 20));
        logMsg('info', msg, 1, 0);
        targetData = filter_hot_pixels(targetData, filterProps{:});
    end

    iTransForm = nTransForms(iFile);
    iRefFrame = nRefFrames(iFile);

    %% reverse
    % in the case of reverse, tform and rframe are the ref -> target
    % tranformation. Therefore, the data/led of the target does not need to
    % be transformed. However, the mask itself needs to be transformed form
    % the reference coordinates to the target coordinates later.

    if kwargs.reverse
        msg = sprintf('reverse doesnt work, yet!');
        logMsg('warn', msg, 1, 0);
        transData = targetData;
        transLed = targetLed;
    else
        msg = sprintf('transforming: target data & LED  << ... %s >>', iFile(end -30:end));
        logMsg('info', msg, 1, 0);
        transData = tform_data(targetData, iTransForm, iRefFrame);
        transLed = tform_data(targetLed, iTransForm, iRefFrame);
    end

    % create struct for the the transformed data of this file
    fileTransForm = struct;

    fileTransForm.refFile = refFile;
    fileTransForm.refData = refData;
    fileTransForm.refLed = refLed;

    fileTransForm.fileName = iFile;
    fileTransForm.targetLed = targetLed;
    fileTransForm.targetData = targetData;

    fileTransForm.transData = transData;
    fileTransForm.transLed = transLed;
    fileTransForm.transForm = iTransForm;
    fileTransForm.refFrame = iRefFrame;
    fileTransForm.reverse = kwargs.reverse;
    fileTransForm.binning = detect_binning(target);

    % save the result in the trans_data container for later use
    transformedData(iFile) = fileTransForm;

    %%% create checkplots
    if kwargs.checkPlot
        if i ~= kwargs.fixedIdx
            check_plot(fileTransForm);
        end
    end
end
end

function check_plot(fileTransForm)
%check_plot(fileTransForm)
figure('units', 'normalized', 'outerposition', [0.2, 0.4, 0.5, 0.5], 'NumberTitle', 'off', 'Name', fileTransForm.fileName);
ax1 = subplot(2, 3, 1);
ref = fileTransForm.refData;
QDM_figure(ref, 'title', 'reference Data', 'ax', ax1)

ax2 = subplot(2, 3, 2);
target = fileTransForm.targetData;
QDM_figure(target, 'title', 'target Data', 'ax', ax2)

ax3 = subplot(2, 3, 3);
trans = fileTransForm.transData;
QDM_figure(trans, 'title', 'transformed Data', 'ax', ax3)

ax4 = subplot(2, 3, 4);
QDM_figure(fileTransForm.refLed, 'title', 'reference LED', 'ax', ax4, 'led', 'true')

ax5 = subplot(2, 3, 5);
QDM_figure(fileTransForm.targetLed, 'title', 'target LED', 'ax', ax5, 'led', 'true')

ax6 = subplot(2, 3, 6);
QDM_figure(fileTransForm.transLed, 'title', 'transformed LED', 'ax', ax6, 'led', 'true')

linkaxes([ax1, ax2, ax3]);
linkaxes([ax4, ax5, ax6]);
end
