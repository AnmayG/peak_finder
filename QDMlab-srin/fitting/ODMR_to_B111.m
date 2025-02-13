function fits = ODMR_to_B111(kwargs)
%[fits] = ODMR_to_B111('nFolders', 'binSizes', 'fieldPolarity', 'type', 'globalFraction', 'forceGuess', 'checkPlot', 'plotGuessSpectra', 'gaussianFit', 'gaussianFilter', 'smoothDegree', 'save', 'savePlots', '['N15',', 'slopeCorrection', 'crop', 'fcrop')
% and then determines B111 field values from the different polarities.
%
% Parameters
% ----------
%   nFolders: cell ['none']
%       if 'none' lets you pick a folder
%   binSize: cell(int) [4]
%   fieldPolarity: (0)
%     0 = neg & pos,
%     1 = neg,
%     2 = pos,
%     4 = nppn
%   globalFraction: numeric (0.25)
%     Amount of global ilumination signal to be subtracted from the
%     measurements before fitting
%   type: [0,1,2] (2)
%     0: for global guess, 
%     1: local with gaussian fit (OLD)
%     2: local with GPU fit NEW
%   smoothDegree: int (2)
%     gaussian smoothing before fit
%   gaussianFit: bool (false)
%   gaussianFilter: numeric (0)
%   nucSpinPol: bool (false)
%   save: bool (true)
%   show diagnostics plots
%   plotGuessSpectra: bool (1)
%   forceGuess: bool (0)
%   slope_correction: bool [false]
%       uses a linear slope correction on the raw data to determine the
%       initial guess. 
%       Note: only works for type = 2
%   crop: bool False
%
%   Returns
%   -------
%   fits: cell (nFolders, binSizes)
%         e.g. fits{1,2} is the first measurement folder with the second binSize
%         e.g. fits{3,1} is the third measurement folder with the first binSize
%
arguments
    kwargs.nFolders {foldersMustExist(kwargs.nFolders)} = 'none';
    kwargs.binSizes = 'none'
    % keyword arguments
    kwargs.fieldPolarity (1,1) {mustBeMember(kwargs.fieldPolarity,[0,1,2,4])} = 0
    kwargs.type (1,1) {mustBeMember(kwargs.type,[0,1,2])} = 2
    kwargs.globalFraction = 'none';
    kwargs.forceGuess (1,1) {mustBeBoolean(kwargs.forceGuess)} = false;
    kwargs.checkPlot (1, 1) {mustBeBoolean(kwargs.checkPlot)} = false;
    kwargs.plotGuessSpectra (1,1) {mustBeBoolean(kwargs.plotGuessSpectra)} = false;
    kwargs.gaussianFit (1,1) {mustBeBoolean(kwargs.gaussianFit)} = false;
    kwargs.gaussianFilter (1,1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0
    kwargs.smoothDegree  (1,1) {mustBeNumeric, mustBePositive} = 2
    
    kwargs.save (1,1) {mustBeBoolean(kwargs.save)} = 1
    kwargs.savePlots (1,1) {mustBeBoolean(kwargs.savePlots)} = 1

    kwargs.diamond {mustBeMember(kwargs.diamond, ...
        ['N15', 'N14', 'DAC', 'singlet', 'doublet', 'triplet', 'gaussian'])} = 'N15'
    kwargs.slopeCorrection = false;
    kwargs.crop  = 'none'
    kwargs.fcrop (1,1) {mustBeBoolean(kwargs.fcrop)}  = false

end
defaults = struct('binSizes', [4], 'globalFraction', 0.25);
nFolders = automatic_input_ui__(kwargs.nFolders);
kwargs = ask_arguments(kwargs, defaults);
binSizes = kwargs.binSizes;

% check if there is more than one folder
nFolders = correct_cell_shape(nFolders);

% check if there is one or more binSize
if isnumeric(binSizes)
    binSizes = [binSizes];
end

% select field polarity
fp = containers.Map({0 1 2 4},{'np  ' 'n   ' 'p   ' 'nppn'});
type = fp(kwargs.fieldPolarity);

fits = cell(size(nFolders,2), size(binSizes, 2));

% set crop to none if 0 or false
if kwargs.crop == 0 | kwargs.crop == false
    kwargs.crop = 'none';
end

for i = 1:size(nFolders,2)
    dataFolder = nFolders{i};
    
    if ~strcmp(kwargs.crop, 'none')
        if kwargs.crop == 1 || kwargs.crop == true
            LED = get_led(dataFolder);
            [~, coordinates] = pick_box(LED, 'led',1, 'title', 'crop Data', 'n',1,'closeFig',true);
            kwargs.crop = coordinates{1};
        end
    end

    for n=1:size(binSizes,2)
        binSize=binSizes(n);
        fit = GPU_fit(dataFolder, binSize,...
                        'fieldPolarity',kwargs.fieldPolarity, ...
                        'type', kwargs.type,...
                        'globalFraction', kwargs.globalFraction,...
                        'diamond', kwargs.diamond,...
                        'gaussianFit', kwargs.gaussianFit,...
                        'gaussianFilter', kwargs.gaussianFilter,...
                        'forceGuess', kwargs.forceGuess,...
                        'checkPlot', kwargs.checkPlot,...
                        'smoothDegree', kwargs.smoothDegree,...
                        'slopeCorrection', kwargs.slopeCorrection,...
                        'crop', kwargs.crop,...
                        'fcrop', kwargs.fcrop,...
                        'save', kwargs.save);
                    
        fit.kwargs = kwargs;
        fit.nFolder = dataFolder;
        
        folderName=[num2str(binSize) 'x' num2str(binSize) 'Binned'];

        fit = get_B111(dataFolder, folderName, type, fit, binSize,...
            'crop', kwargs.crop);
        
        if kwargs.savePlots
            ODMR_to_B111_plot(fit, fullfile(dataFolder, folderName));
        end

        if kwargs.save
            % copy laser image and csv
            msg = sprintf('copying laser.jpg, laser.csv into %s', dataFolder);
            logMsg('info',msg,1,0);

            copyfile(fullfile(dataFolder, 'laser.csv'),fullfile(dataFolder, folderName), 'f')
            copyfile(fullfile(dataFolder, 'laser.jpg'),fullfile(dataFolder, folderName), 'f')

            fName = sprintf('final_fits_(%ix%i).mat', binSize, binSize);
            msg = sprintf('saving %s into %s', fName, dataFolder);
            logMsg('info',msg,1,0);

            save(fullfile(dataFolder, folderName, fName), '-struct', 'fit');
        end
        fits{i,n} = fit;
    end
end

end
