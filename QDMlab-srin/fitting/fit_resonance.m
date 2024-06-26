function fit = fit_resonance(expData, kwargs)
    %[fit] = fit_resonance(expData, binSize; 'type', 'checkPlot', 'gaussianFilter', 'model', 'slopeCorrection', 'crop', 'fcrop')
    % fits a single resonance frequency (i.e. low/high frequency range) of
    % either positive or negative field.
    %
    % ----
    %  This is what the function does:
    %
    %  1. Extract signal and freq from expData
    %  2. global_guess if `type` ~= 2
    %  3. reshapes data for gpufit -> all pixels in a row
    %  4. if *type* == 1: guess peaks; if *type* == 2: **(a)** get_initial_guess -> creates pre guess for single gaussian GPU fit; **(b)** gpu_fit (GAUSS_1D); **(c)** parameters_to_guess calculates the initial guess from the fitted parameters of 6.b.
    %  5. :code:`gpu_fit` calculates lorentzian fits
    %  6. :code:`reshape_fits` creates the (y,x) sized array out of the fitted parameters
    %  7. *checkPlot* of fits -> needs to be closed to proceed
    %
    % Parameters
    % ----------
    %     expData: struct
    %         Data of :code:`load(run0000n.mat)`
    %     type: int (2)
    %         type of initial guess:
    %         0: global
    %         1: local
    %         2: gaussian
    %     checkPlot: int (0)
    %         Creates an interactive plot to check the fits
    %  state definitions: CONVERGED = 0, MAX_ITERATION = 1,
    %                     SINGULAR_HESSIAN = 2, NEG_CURVATURE_MLE = 3,
    %                     GPU_NOT_READY = 4
    
    arguments
        expData struct
        kwargs.type (1, 1) {mustBeMember(kwargs.type, [0, 1, 2])} = 0; %%2
        kwargs.checkPlot (1, 1) {mustBeBoolean(kwargs.checkPlot)} = false;
        kwargs.model {mustBeMember(kwargs.model, ['N15', 'N14', 'DAC', 'singlet', 'doublet', 'gaussian', 'quad', 'gauss4'])} = 'gauss4';%'doublet';
        kwargs.slopeCorrection = false;
    end
    show_references()
    
    msg = sprintf('--------------------------------------------------------------------');
    logMsg('info',msg,1,0);
    tStart = tic;
    
    % preallocate fit
    fit = struct();
    
    %% data preparation
    binData = expData.signal;
    freq = expData.SweepParam / 1e9;
    sizeX = size(binData, 2); % binned image x-dimensions
    sizeY = size(binData, 1); % binned image y-dimensions
    binDataNorm = zeros(size(binData));
    NormalizationFactor = mean(binData,3); % get average of intensities
    for y = 1:length(freq)
        binDataNorm(:,:,y) = binData(:,:,y) ./ NormalizationFactor;
    end

    msg = sprintf('starting parameter estimation (%s)', kwargs.model');
    logMsg('info',msg,1,0);
    
    %% Prepare gpuData for fitting in GPU
    sweepLength = size(freq, 2); % srin - i think this is correct because its 1xn regardless
    imgPts = sizeX * sizeY; % number of (x,y) pixels
    gpudata = reshape(binDataNorm, [imgPts, sweepLength]); % make it into 2d matrix
    gpudata = transpose(gpudata);
    gpudata = single(gpudata);
    if kwargs.slopeCorrection
        gpudata = slope_correction(gpudata, freq, kwargs.slopeCorrection);
    end
    xValues = single(freq');

    %% GUESS INITIAL FIT PARAMETERS
    
    %% Determine global guess if not given one
    if kwargs.type ~= 2
        initialGuess = global_guess(binDataNorm, freq); % initial guess for GPUfit
    end
    
    %% global guess
    if kwargs.type == 0 % reshape into [6 x numpoints]
        initialGuess = reshape(initialGuess, [imgPts, 13]);
        initialGuess = transpose(initialGuess);
    end
    
    %% local guess -> guess parameter for each pixel
    
    %% output setup
    fit.initialGuess.states = zeros(size(binDataNorm, [1 2]));
    if kwargs.type == 1 %% old local/gaussian guess
        msg = sprintf('local guess estimation');
        logMsg('info',msg,1,0);
        sizeX = size(binDataNorm, 2); % binned image x-dimensions
        sizeY = size(binDataNorm, 1); % binned image y-dimensions
        freqstep=freq(2)-freq(1);

        %% generating the initialguess for fitting
        % iterate over all pixels with index (x,y)
        h = waitbar(0, "Starting local guess");
        split = 3.1e9;
        tolerancepeak = 1e8;
        pixelData = binDataNorm;
        freq2=freq*1E9; %peakfinder currently uses Hz, GPU uses GHz
        parfor x = 1:sizeX
            % waitbar(x / sizeX, h, "Column " + x);
            for y = 1:sizeY
                seed=peak_find_function(false, split, tolerancepeak, y, x, pixelData, freq2);
                initialGuess(y, x, :) = seed(1:13);
            end
        end
        close(h)
        initialGuess = reshape(initialGuess, [imgPts, 13]);
        initialGuess = transpose(initialGuess);
    end
    
    %% GPU pre fits
    
    %% fitting related
    tolerance = 1e-13;
    initialPreGuess = 'none';
    
    if kwargs.type == 2
        %% initial preGuess
        initialPreGuess = get_initial_guess(gpudata, freq, kwargs.model);
        fit.initialPreGuess = initialPreGuess;
        
        % N15 (i.e. doublet) use the preGuess instead of a gaussian
        if find(strcmp(kwargs.model, {'N15', 'doublet'}))
            msg = sprintf('determining initial guess only from (N14) preInitialGuess'); 
            logMsg('debug',msg,1,0);
            
            % Note: DAC initial guess is different than N15
            initialGuess = parameters_to_guess(initialPreGuess, kwargs.model);
            fit.initialGuess.states = zeros(size(initialGuess));
        else
            % single gaus fit for initial parameters
            model_id = ModelID.GAUSS_1D;
            [initialGuess, states, chiSquares, ~, ~] = gpufit(gpudata, [], ...
                model_id, initialPreGuess, tolerance, 1000, ...
                [], EstimatorID.MLE, xValues);
            initialGuess = parameters_to_guess(initialGuess, kwargs.model);
            fit.initialGuess.chi = chiSquares;
            fit.initialGuess.states = states;
        end
    end
    
    nBadPixels = numel(nonzeros(fit.initialGuess.states));
    
    if nBadPixels > 1
        if kwargs.type == 2 && strcmp(kwargs.model, 'N14')
            msg = sprintf('%i / %i pixels failed the pre-guess. See fits.states', nBadPixels, imgPts);
            logMsg('warn',msg,1,0);
        else
            msg = sprintf('%i / %i pixels had to be substituted', nBadPixels, imgPts);
            logMsg('warn',msg,1,0);
        end
    end
    msg = sprintf('%i: initial parameter estimation complete in: %.1f s', toc(tStart));
    logMsg('info',msg,1,0);
    
    %% FINAL GPU FIT
    switch kwargs.model
        case {'gauss4'}
            model_id = ModelID.GAUSS4;
        case {'N14', 'triplet'}
            model_id = ModelID.ESR14N;
        case {'N15', 'doublet'}
            model_id = ModelID.ESR15N;
        case {'DAC','singlet'}
            model_id = ModelID.ESRSINGLE;
        case 'gaussian'
            model_id = ModelID.GAUSS_1D;
        case 'lin'
            model_id=ModelID.LINEAR_1D;
        case 'quad'
            model_id=ModelID.QUAD_2;
    end
    
    max_n_iterations = 500;
    msg = sprintf('starting GPU fit, model: %s', kwargs.model);
    logMsg('info',msg,1,0);
    
    %BELOW, the standard pattern for low pressure should be stds, devs, stds fitting
    %REPLACE gpuweight with [] if you dont want to use weights, or vice versa,
    %it is the second parameter in gpufit (only on option 1 currently)
    [parameters, ~, ~, ~, ~] = gpufit(gpudata, [], ...
        model_id, initialGuess, tolerance, max_n_iterations, [1,1,1,0,1,0,0,0,1,0,1,0,1].', EstimatorID.LSE, xValues);
    msg = sprintf('Round 1 Done');
    logMsg('info',msg,1,0);
    
    [parameters, states, chiSquares, n_iterations, ~] = gpufit(gpudata, [], ...
        model_id, parameters, tolerance, max_n_iterations, [0,0,0,1,0,1,1,1,0,1,0,1,0].', EstimatorID.LSE, xValues);
    
    
    
    % failed fits for pixel with extrem chiSquared or values outside of the
    % measurement freq. range
    % states(chiSquares > 5e-4) = 5;
    
    fit = make_fit_struct(fit, initialPreGuess, initialGuess, parameters, states, ...
        chiSquares, n_iterations, sizeX, sizeY, kwargs.model);
    fit.freq = freq;
    
    msg = sprintf('final GPU fitting complete in: %.1f s', toc(tStart));
    logMsg('info',msg,1,0);
    
    if numel(nonzeros(states)) > 0
        badPre = numel(nonzeros(states));
        msg = sprintf('%i / %i pixels failed the final fit. See fit.states!', badPre, imgPts);
        logMsg('warn',msg,1,0);
    end
    
    if kwargs.checkPlot
        msg = sprintf('close figure to continue');
        logMsg('ATTENTION',msg,1,0);
        fig = gpu_fit_checkPlot(fit, binDataNorm, freq, 1, kwargs.model);
        waitfor(fig)
    end
end

function data = slope_correction(data, freq, nPoints)
    % [data] = slope_correction(data, freq, nPoints)
    % calculates slope between 1st - last pixel and removes it from data
    % for easy fitting
    msg = sprintf('correcting slope of for the initial guess calculation');
    logMsg('debug',msg,1,0);
    
    d1 = nanmean(data(1:nPoints+1,:));
    dend = nanmean(data(end+1-nPoints:end,:));
    
    delta = (dend-d1);
    slope = delta/numel(freq);
    
    correction = zeros(size(data));

    for i = 1:numel(freq)
        correction(i,:) = slope*(i-1);
    end

    data = data - correction;
end

%% fitting helper functions
function initialGuess = get_initial_guess(gpudata, freq, diamond)
    initialGuess = zeros(4, size(gpudata, 2), 'single');
    
    % amplitude
    mx = max(gpudata, [], 'omitnan');
    mn = min(gpudata, [], 'omitnan');
    initialGuess(1, :) = -abs(((mx - mn)./mx));
    
    % center frequency
    [~, idx] = sort(gpudata);
    l = 10; % lowest n values
    mxidx = max(idx(1:l, :));
    mnidx = min(idx(1:l, :));
    
    switch diamond
        case {'N15', 'doublet', 'gaussian', 'DAC', 'singlet', 'gauss4'}
            cIdx = int16((mxidx+mnidx)/2);
        case {'N14', 'triplet'}
            cIdx = int16(mean(cat(1, mxidx, mnidx), "omitnan"));
    end
    
    center = freq(cIdx);
    initialGuess(2, :) = center;
    
    % width
    if strcmp(diamond, 'DAC')
        initialGuess(3, :) = 0.004;
    elseif strcmp(diamond, 'singlet')
        initialGuess(3, :) = 0.0002;
    else
        initialGuess(3, :) = 0.0005;
    end
        
    % offset
    initialGuess(4, :) = mx;
    end
    
    function guess = parameters_to_guess(parameters, diamond)
    %[guess] = parameters_to_guess(parameters, diamond)
    guess(1, :) = parameters(2, :); % location
    guess(2, :) = 0.0002; % width
    switch diamond
        case {'N14', 'triplet'}
            guess(3, :) = abs(parameters(1, :)); % amplitude (contrast)
            guess(4, :) = abs(parameters(1, :)); % amplitude (contrast)
            guess(5, :) = abs(parameters(1, :)); % amplitude (contrast)
            guess(6, :) = parameters(4, :) - 1; % baseline
        case {'N15', 'doublet'}
            guess(2, :) = 0.0004; % width
            guess(3, :) = -parameters(1, :); % amplitude (contrast)
            guess(4, :) = -parameters(1, :); % amplitude (contrast)
            guess(5, :) = parameters(4, :) - 1; % baseline
        case {'singlet'}
            guess(2, :) = 0.0003; % width
            guess(3, :) = -parameters(1, :); % amplitude (contrast)
            guess(4, :) = parameters(4, :)-1; % baseline
        case {'DAC'}
            guess(2, :) = 0.005; % width
            guess(3, :) = -parameters(1, :); % amplitude (contrast)
            guess(4, :) = parameters(4, :)-1; % baseline
        case 'gaussian'
            guess = parameters;
       guess = single(guess);
    end
end

function fit = make_fit_struct(fit, preGuess, initialGuess, parameters, states, chiSquares, n_iterations, sizeX, sizeY, diamond)
    msg = sprintf('reshaping data into (%4i, %4i)', sizeY, sizeX);
    logMsg('debug',msg,1,0);
    
    %make parameters matrix into 3d matrix with x pixels, y pixels, and parameters
    %     fit.preGuess = reshape(preGuess, [], sizeY, sizeX); % initial guess
    fit.initialGuess.parameters = double(reshape(initialGuess, [], sizeY, sizeX)); % initial guess
    fit.initialGuess.p = double(initialGuess);
    
    fit.parameters = double(reshape(parameters, [], sizeY, sizeX)); % fitted parameters
    fit.p = double(parameters);
    
    % matricies with 2 dimensions for x and y pixels:
    fit.resonance = squeeze(fit.parameters(1, :, :));
    fit.width = squeeze(fit.parameters(2, :, :));
    %%fit.contrastA = squeeze(fit.parameters(3, :, :));
    
    % check for diamond type: N14 has 6 parameters, N15 has only 5
    nParams = size(initialGuess, 1);
    
    switch diamond
        case {'N14', 'triplet'} % for N14 diamonds
            fit.contrastB = squeeze(fit.parameters(4, :, :));
            fit.contrastC = squeeze(fit.parameters(5, :, :));
            fit.baseline = squeeze(fit.parameters(6, :, :)+1);
        case {'N15', 'doublet'} % for N15 diamonds
            fit.contrastB = squeeze(fit.parameters(4, :, :));
            fit.baseline = squeeze(fit.parameters(5, :, :)+1);
        case {'singlet', 'DAC'} % singlet
            fit.baseline = squeeze(fit.parameters(4, :, :)+1);
        case {'gaussian'}
            fit.resonance = squeeze(fit.parameters(2, :, :));
            fit.width = squeeze(fit.parameters(3, :, :));
            fit.contrastA = -squeeze(fit.parameters(1, :, :));
            fit.baseline = squeeze(fit.parameters(4, :, :));
    end
    
    fit.states = reshape(states, [sizeY, sizeX]);
    fit.chiSquares = double(reshape(chiSquares, [sizeY, sizeX]));
    fit.n_iterations = reshape(n_iterations, [sizeY, sizeX]);
    fit.binSize = 1;
    
    if ~strcmp(preGuess, 'none')
        fit.pg = double(parameters_to_guess(preGuess, diamond));
    else
        fit.pg = preGuess;
    end
end
