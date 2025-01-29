function fit = fit_resonance_pf(sig, sp, pdf, binSize, nRes, header, kwargs)
%[fit] = fit_resonance(expData, binSize, nRes; 'header', 'type', 'globalFraction', 'forceGuess', 'checkPlot', 'gaussianFit', 'gaussianFilter', 'smoothDegree', 'diamond', 'slopeCorrection', 'crop', 'fcrop')
% fits a single resonance frequency (i.e. low/high frequency range) of
% either positive or negative field.
%
% hint
% ----
%  This is what the function does:
%
%  1. :code:`prepare_raw_data`
%  2. gaussian filter (:code:`imfilter`) if `gaussianFilter` == 1
%  3. Global illumination correction with :code:`correct_global`
%  4. global_guess if `type` ~= 2
%  5. reshapes data for gpufit -> all pixels in a row
%  6. if *type* == 1: guess peaks; if *type* == 2: **(a)** get_initial_guess -> creates pre guess for single gaussian GPU fit; **(b)** gpu_fit (GAUSS_1D); **(c)** parameters_to_guess calculates the initial guess from the fitted parameters of 6.b.
%  7. :code:`gpu_fit` calculates lorentzian fits
%  8. :code:`reshape_fits` creates the (y,x) sized array out of the fitted parameters
%  9. *checkPlot* of fits -> needs to be closed to proceed
%
%
% Parameters
% ----------
%     expData: struct
%         Data of :code:`load(run0000n.mat)`
%     binSize: int
%         binning size (can be 1)
%     nRes: int
%         number of resonance. Low frequencies = 1, High frequencies = 2
%     type: int (2)
%         type of initial guess:
%         1: global
%         2: local
%         3: gaussian
%     globalFraction: double (0.5)
%         Ammount of global signal to be corrected for (see. correct_global)
%     forceGuess: int (0)
%         Used for forcing a guess (NOT IMPLEMENTED)
%     checkPlot: int (0)
%         Creates an interactive plot to check the fits
%     gaussianFit: int (0)
%         In case the :code:`type = local` and the MATLAB function find_peaks does
%         not find 3 peaks.
%         **if 0**: the global guess will be used for that pixel
%         **if 1**: a gaussian fit is used to find peak positions
%     gaussianFilter: int (0)
%         Determines if a gaussian filter is applied before fitting
%     smoothDegree: int (2)
%         The ammount of smoothing if gaussianFilter == 1
%     diamond: str (N14)
%         The type of diamond. Choses the type of fitting.
%
%  state definitions: CONVERGED = 0, MAX_ITERATION = 1,
%                     SINGULAR_HESSIAN = 2, NEG_CURVATURE_MLE = 3,
%                     GPU_NOT_READY = 4,

arguments
    sig
    sp
    pdf
    binSize double
    nRes (1, 1) int16
    header = 'none'
    kwargs.type (1, 1) {mustBeMember(kwargs.type, [0, 1, 2])} = 0; %%2
    kwargs.globalFraction (1, 1) {mustBeNumeric} = 0 %0.5 formerly srin
    kwargs.forceGuess (1, 1) {mustBeMember(kwargs.forceGuess, [1, 0])} = false;
    kwargs.checkPlot (1, 1) {mustBeBoolean(kwargs.checkPlot)} = false;
    kwargs.gaussianFit (1, 1) {mustBeBoolean(kwargs.gaussianFit)} = false;
    kwargs.gaussianFilter (1, 1) {mustBeNumeric, mustBeGreaterThanOrEqual(kwargs.gaussianFilter, 0)} = 0;
    kwargs.smoothDegree (1, 1) {mustBeNumeric, mustBePositive} = 2
    kwargs.diamond {mustBeMember(kwargs.diamond, ['N15', 'N14', 'DAC', 'singlet', 'doublet', 'gaussian', 'quad', 'gauss4', 'gauss2', 'gaussTwin', 'gaussTwinD'])} = 'gaussTwinD'; %'gauss4';%'doublet';
    kwargs.slopeCorrection = false;
    kwargs.crop = 'none'
    kwargs.fcrop  = false

end
show_references()

msg = sprintf('--------------------------------------------------------------------');
logMsg('info',msg,1,0);
tStart = tic;

% preallocate fit
fit = struct();

%% check type/diamond combination
if kwargs.type ~= 2 && ~any(strcmpi(kwargs.diamond, {'N14', 'triplet'}))
    msg = sprintf('Determining the initial parameters for a fit with this method is not supported for N15 diamonds, yet');
    logMsg('error',msg,1,0);
end

%% data preparation
% this step could easily be skipped, the only thing one needs to figure out
% is how to get the
[binDataNorm, freq] = prepare_raw_data_pf(sig,sp, binSize, nRes, header, 'crop', kwargs.crop, 'fcrop', kwargs.fcrop);

%%srin TEST an array of all 1.07

%binDataNorm=ones(1, 1, 400) * -7.6485;

sizeX = size(binDataNorm, 2); % binned image x-dimensions
sizeY = size(binDataNorm, 1); % binned image y-dimensions

%% 1. GAUSSIAN BLUR
% default = 0
% gaussian filter on the guesses, can lead to some problems if high gradient
%%SRIN COMMENTED THIS OUT, currently blurring within peakfind
% if kwargs.gaussianFilter ~= 0
%     msg = sprintf('%i: smoothing data using gaussian blur: %.1f', nRes, gaussianFilter');
%     logMsg('info',msg,1,0);
%     gFilter = fspecial('gaussian', [20, 20], gaussianFilter);
%     binDataNorm = imfilter(binDataNorm, gFilter, 'symmetric', 'conv');
% end
%%END srin
msg = sprintf('%i: starting parameter estimation (%s)', nRes, kwargs.diamond');
logMsg('info',msg,1,0);

%% global spectra subtraction

%%srin currently correct_global is BREAKING THE DATA and making it all NaN!
%binDataNorm = correct_global(binDataNorm, kwargs.globalFraction); %%SRIN
%commented this

%% first determine global guess
meanData = squeeze(mean(binDataNorm, [1, 2], 'omitnan'));

% if kwargs.type ~= 2
%     initialGuess = global_guess(binDataNorm, freq); % initial guess for GPUfit
% end

%% prepare GPUfit data
sweepLength = size(freq, 2); %%srin i think this is correct because its 1xn regardless
imgPts = sizeX * sizeY; % number of (x,y) pixels

gpudata = reshape(binDataNorm, [imgPts, sweepLength]); % make it into 2d matrix
gpudata = transpose(gpudata); %transpose to make it 51 x pixels

gpudata = single(gpudata);

if kwargs.slopeCorrection
    gpudata = slope_correction(gpudata, freq, kwargs.slopeCorrection);
end

xValues = single(freq');
%disp(xValues)
%% GUESS INITIAL FIT PARAMETERS

%% gloabl guess
% if kwargs.type == 0 % reshape into [6 x numpoints]
%     initialGuess = reshape(initialGuess, [imgPts, 13]); %%srinsrin originally 13
%     initialGuess = transpose(initialGuess);
% end

%% local guess -> guess parameter for each pixel

%% output setup
fit.initialGuess.states = zeros(size(binDataNorm, [1 2]));

if kwargs.type == 1 %% old local/gaussian guess
    msg = sprintf('%i: local guess estimation', nRes');
    logMsg('info',msg,1,0);
    
    
    sizeX = size(binDataNorm, 2); % binned image x-dimensions
    sizeY = size(binDataNorm, 1); % binned image y-dimensions
    weights=ones(sizeY,sizeX,sweepLength);
    freqstep=freq(2)-freq(1);
    weightbin=0.02;
    weightrange=int32(weightbin/freqstep);

    %% generating the initialguess for fitting
    % iterate over all pixels with index (x,y)
    for x = 1:sizeX
        for y = 1:sizeY %eventually, remake peakfinding with gpu arrays
            split = 3.1e9;
            tolerancepeak = 0.1e9;
            %pixelData = squeeze(binDataNorm(y, x, :));
            pixelData = binDataNorm;
            tempweight=zeros(sweepLength,1);
            %disp(size(pixelData));
            freq2=freq*1E9; %peakfinder currently uses Hz, GPU uses GHz
            seed=peak_find_function(split, tolerancepeak, y, x, pixelData, freq2);
            %find(abs(freq-2.7617)<0.001)
            
            fndpks=[find(abs(freq-(seed(1)-seed(2)))<0.001), find(abs(freq-(seed(1)+seed(2)))<0.001)]; %lol just temporary just return the peak index
            for idx = fndpks
                tempweight(max(1, idx-weightrange):min(length(tempweight), idx+weightrange)) = 1; %so only fit data within 50 MHz of each found peak
            end
            



            if abs(seed(1)-split/1E9)<5*tolerancepeak/1E9 %i.e. if your fit D value is reasonable, I take the findpeaks to be 'good' and update the local seed, else, default to keeping the global_guess seed. the 5*tolerance threshold is arbitrary and there should be an interal flag for good peak finding in the future. 
                initialGuess(y, x, :) = seed;
                
                weights(y,x,:) = tempweight;
            end

            % if x==91 & y==91
            %     disp(weights(y,x,:))
            % end

            %% SRIN, old option 1 code from Fu

            % [pkVal, pkLoc, fitFlg] = guess_peaks(pixelData, meanData, freq, ...
            %     'smoothDegree', kwargs.smoothDegree, ...
            %     'forceGuess', kwargs.forceGuess, ...
            %     'gaussianFit', kwargs.gaussianFit, ...
            %     'pixel', [y, x, nRes]);
            % 
            % % check if find peaks returned 3 peaks
            % % add them to fit.initialGuess.states
            % if fitFlg ~= 0 % 1 == gauss 2= global (i.e. local failed)
            %     fit.initialGuess.states(y,x) = fitFlg;
            % end
            % 
            % % replace guess with new guess if fitFlg is not global
            % % (i.e. fitFlg = 2)
            % if fitFlg ~= 2
            %     % if it returns 3 -> replace the global guess with
            %     % the local
            %     resonance = (pkLoc(1) + pkLoc(2) + pkLoc(3)) / 3; % in GHz
            %     width = 0.0005;
            %     contrast = (nanmean(pixelData(1:10)) + pkVal - 1)';
            %     baseline = nanmean(pixelData(1:10)) - 1;
            %     initialGuess(y, x, :) = [resonance, width, contrast, baseline];
            % end
            %%end SRIN
        end
    end
    initialGuess = reshape(initialGuess, [imgPts, 13]);
    initialGuess = transpose(initialGuess);

    gpuweight = reshape(weights, [imgPts, sweepLength]); % make it into 2d matrix
    gpuweight = transpose(gpuweight); 
    gpuweight = single(gpuweight);    

end

%% GPU pre fits

%% fitting related
tolerance = 1e-13;

initialPreGuess = 'none';

if kwargs.type == 2
    %% initial preGuess
    initialPreGuess = get_initial_guess(gpudata, freq, kwargs.diamond);
    fit.initialPreGuess = initialPreGuess;
    
    % N15 (i.e. doublet) use the preGuess instead of a gaussian
    if find(strcmp(kwargs.diamond, {'N15', 'doublet'}))
        msg = sprintf('determining initial guess only from (N14) preInitialGuess'); 
        logMsg('debug',msg,1,0);
        
        % Note: DAC initial guess is different than N15
        initialGuess = parameters_to_guess(initialPreGuess, kwargs.diamond);
        fit.initialGuess.states = zeros(size(initialGuess));
    else
        % single gaus fit for initial parameters
        model_id = ModelID.GAUSS_1D;
        [initialGuess, states, chiSquares, n_iterations, time] = gpufit(gpudata, [], ...
            model_id, initialPreGuess, tolerance, 1000, ...
            [], EstimatorID.MLE, xValues);
        initialGuess = parameters_to_guess(initialGuess, kwargs.diamond);
        fit.initialGuess.chi = chiSquares;
        fit.initialGuess.states = states;
    end
end

nBadPixels = numel(nonzeros(fit.initialGuess.states));

if nBadPixels > 1
    if kwargs.type == 2 && strcmp(kwargs.diamond, 'N14')
        msg = sprintf('%i: %i / %i pixels failed the pre-guess. See fits.states', nRes, nBadPixels, imgPts);
        logMsg('warn',msg,1,0);
    else
        msg = sprintf('%i: %i / %i pixels had to be substituted', nRes, nBadPixels, imgPts);
        logMsg('warn',msg,1,0);
    end
end
msg = sprintf('%i: initial parameter estimation complete in: %.1f s', nRes, toc(tStart));
logMsg('info',msg,1,0);

%% FINAL GPU FIT
switch kwargs.diamond
    case {'gauss4'}
     model_id = ModelID.GAUSS4;
    case {'gauss2'}
     model_id = ModelID.GAUSS2;
    case {'gaussTwin'}
     model_id = ModelID.GAUSSTWIN;
    case {'gaussTwinD'}
     model_id = ModelID.GAUSSTWIND; 
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

%%srin the max iterations was originally 1000, I updated it to 100000
max_n_iterations = 1000;
%disp(gpudata);
msg = sprintf('%i: starting GPU fit, model: %s', nRes);
logMsg('info',msg,1,0);
%[1,1,1,1,1,1,1,1,1,1,1,1,1].'
%[1,1,1,0,1,0,0,0,0,0,0,0,1].' 
%[1,1,1,0,1,0,0,0,1,0,1,0,1].' pattern of stds
%[0,0,0,1,0,1,1,1,0,1,0,1,0].' pattern of devs
%[0,0,0,0,0,0,0,0,0,0,0,0,1].'
% run Gpufit
%%srin previously EstimatorID.MLE


%%SRIN COMMENT BLOCK: GPUFIT WITH CONSTRAINT
% 
% constraint = zeros(sizeY, sizeX, 26);
% constraint(:,:,1) = 2.85;
% constraint(:,:,2) = 3;
% 
% constraint(:,:,3) = 0;
% constraint(:,:,4) = 0.4;
% 
% constraint(:,:,5) = 0;
% constraint(:,:,6) = 0.05;
% 
% constraint(:,:,7) = -0.015;
% constraint(:,:,8) = 0.015;
% 
% constraint(:,:,9) = -0.04;
% constraint(:,:,10) = 0;
% 
% constraint(:,:,11) = -0.01;
% constraint(:,:,12) = 0.01;
% 
% constraint(:,:,13) = 0;
% constraint(:,:,14) = 0.01;
% 
% constraint(:,:,15) = -0.01;
% constraint(:,:,16) = 0.01;
% 
% constraint(:,:,17) = 0;
% constraint(:,:,18) = 0.05;
% 
% constraint(:,:,19) = -0.012;
% constraint(:,:,20) = 0.012;
% 
% constraint(:,:,21) = -0.03;
% constraint(:,:,22) = 0;
% 
% constraint(:,:,23) = -0.01;
% constraint(:,:,24) = 0.01;
% 
% constraint(:,:,25) = 0.98;
% constraint(:,:,26) = 1.02;
% 
% constraint = reshape(constraint, [imgPts, 26]); %%srinsrin originally 13
% initconstraint = single(transpose(constraint));
% contype= int32(ones(1, 13) * 3); %%srin i.e. constrain both upper and lower bounds, see constrainttype.m
% 
% [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
%     model_id, initialGuess, initconstraint, contype, tolerance, max_n_iterations, [1,1,1,0,1,0,0,0,1,0,1,0,1].', EstimatorID.LSE, xValues);
% disp("round 1 done")
% 
% [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
%     model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [0,0,0,1,0,1,1,1,0,1,0,1,0].', EstimatorID.LSE, xValues);

%%END SRIN COMMENT BLOCK: GPUFIT WITH CONSTRAINT


% %%SRIN COMMMENT BLOCK: GPUFIT WITHOUT CONSTRAINT



%%srin now run it again using the fits as initial guesses and allowing
%%everything to vary
%%srin make sure to try implementing gpufit_constrained. the way
%%constraints work is that you give a lower and upper for each one, but by
%%choosing constraint type you get to pick which (neither, upper, lower, or
%%both) are enforced


%high p [1,1,1,0,1,0,0,0,0,0,0,0,1].'
%USE FOR GAUSS4
% [parameters, states, chiSquares, n_iterations, ~] = gpufit(gpudata, [], ...
%     model_id, initialGuess, tolerance, max_n_iterations, [1,1,0,0,0,0,0,0,0,0,0,0,1].', EstimatorID.LSE, xValues);
% disp('round 1')
% [parameters, states, chiSquares, n_iterations, ~] = gpufit(gpudata, [], ...
%     model_id, parameters, tolerance, max_n_iterations, [0,0,1,1,0,0,0,0,0,0,0,0,1].', EstimatorID.LSE, xValues);
% disp ('round 2')
% [parameters, states, chiSquares, n_iterations, ~] = gpufit(gpudata, [], ...
%     model_id, parameters, tolerance, max_n_iterations, [0,0,0,0,1,1,0,0,0,0,0,0,1].', EstimatorID.LSE, xValues);


% [parameters, states, chiSquares, n_iterations, ~] = gpufit(gpudata, [], ...
%     model_id, initialGuess, tolerance, max_n_iterations, [1,1,0,0,0,0,0,0,0,0,0,1].', EstimatorID.LSE, xValues);
% disp('round 1')
% [parameters, states, chiSquares, n_iterations, ~] = gpufit(gpudata, [], ...
%     model_id, parameters, tolerance, max_n_iterations, [0,0,1,1,0,1,0,1,0,1,0,1].', EstimatorID.LSE, xValues);
% disp('round 2')
% [parameters, states, chiSquares, n_iterations, ~] = gpufit(gpudata, [], ...
%     model_id, parameters, tolerance, max_n_iterations, [0,0,1,1,1,1,1,1,1,1,1,0].', EstimatorID.LSE, xValues);

% disp ('round 2')
% [parameters, states, chiSquares, n_iterations, ~] = gpufit(gpudata, [], ...
%     model_id, parameters, tolerance, max_n_iterations, [0,0,0,0,1,1,0,0,0,0,0,0,1].', EstimatorID.LSE, xValues);



%BELOW, the standard pattern for low pressure should be stds, devs, stds fitting
%REPLACE gpuweight with [] if you dont want to use weights, or vice versa,
%it is the second parameter in gpufit (only on option 1 currently)

% [parameters, states, chiSquares, n_iterations, ~] = gpufit(gpudata, gpuweight, ...
%     model_id, initialGuess, tolerance, max_n_iterations, [1,1,1,0,1,0,0,0,1,0,1,0,1].', EstimatorID.LSE, xValues);
% disp("round 1 done")
% 
% [parameters, states, chiSquares, n_iterations, ~] = gpufit(gpudata, gpuweight, ...
%     model_id, parameters, tolerance, max_n_iterations, [0,0,0,1,0,1,1,1,0,1,0,1,0].', EstimatorID.LSE, xValues);

% [parameters, states, chiSquares, n_iterations, ~] = gpufit(gpudata, [], ...
%     model_id, parameters, tolerance, max_n_iterations, [1,1,1,1,1,1,1,1,1,1,1,1,1].', EstimatorID.LSE, xValues);



% %%END SRIN COMMMENT BLOCK: GPUFIT WITHOUT CONSTRAINT


% %%BEGIN SRIN COMMENT BLOCK: GAUSSTWIN WITH CONSTRAINT
% 
% constraint = zeros(sizeY, sizeX, 24);
% constraint(:,:,1) = 3;
% constraint(:,:,2) = 3.2;
% 
% constraint(:,:,3) = 2.8e-3*150;
% constraint(:,:,4) = 2.8e-3*250;
% 
% constraint(:,:,5) = 0;
% constraint(:,:,6) = 0.1;
% 
% constraint(:,:,7) = -0.03;
% constraint(:,:,8) = -0.0005;
% 
% constraint(:,:,9) = -0.0015; %for g1 the height gaps can be up to 0.015
% constraint(:,:,10) = 0.0015;
% 
% constraint(:,:,11) = -0.03;
% constraint(:,:,12) = -0.0005;
% 
% constraint(:,:,13) = -0.0015; %made it rlly small (prev 0003
% constraint(:,:,14) = 0.0015;
% 
% constraint(:,:,15) = 0.005;
% constraint(:,:,16) = 0.015;
% 
% constraint(:,:,17) = -0.005;
% constraint(:,:,18) = 0.005;
% 
% constraint(:,:,19) = 0.005;
% constraint(:,:,20) = 0.015;
% 
% constraint(:,:,21) = -0.005;
% constraint(:,:,22) = 0.005;
% 
% constraint(:,:,23) = 0.98;
% constraint(:,:,24) = 1.02;
% 
% 
% 
% constraint = reshape(constraint, [imgPts, 24]); %%srinsrin originally 13
% initconstraint = single(transpose(constraint));
% contype= int32(ones(1, 12) * 3); %%srin i.e. constrain both upper and lower bounds, see constrainttype.m
% 
% [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
%     model_id, initialGuess, initconstraint, contype, tolerance, max_n_iterations, [1,1,0,0,0,0,0,0,0,0,0,0].', EstimatorID.LSE, xValues);
% disp('round 1')
% [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
%     model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [0,0,1,1,0,1,0,1,0,1,0,0].', EstimatorID.LSE, xValues);
% disp('round 2')
% [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
%     model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [0,0,0,0,1,0,1,0,1,0,1,0].', EstimatorID.LSE, xValues);
% % disp('round 3')
% % [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
% %     model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [1,1,1,0,1,0,1,0,1,0,1,0].', EstimatorID.LSE, xValues);
% % disp('round 4')
% % [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
% %     model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [0,0,1,1,0,1,0,1,0,1,0,0].', EstimatorID.LSE, xValues);
% disp('round 5')
% [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
%     model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [1,1,1,0,0,0,0,0,0,0,0,0].', EstimatorID.LSE, xValues);
% disp('round 5')
% [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
%     model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [1,1,1,1,1,1,1,1,1,1,1,0].', EstimatorID.LSE, xValues);
% 

%%BEGIN SRIN COMMENT BLOCK: GAUSSTWIN WITH CONSTRAINT

constraint = zeros(sizeY, sizeX, 26);
constraint(:,:,1) = 3;
constraint(:,:,2) = 3.1;

constraint(:,:,3) = 2.8e-3*190;
constraint(:,:,4) = 2.8e-3*210;

constraint(:,:,5) = -0.02;
constraint(:,:,6) = 0.02; %prev 0.05 span

constraint(:,:,7) = -0.03;
constraint(:,:,8) = -0.0001;

constraint(:,:,9) = -0.03; %for g1 the height gaps can be up to 0.015
constraint(:,:,10) = -0.0001;

constraint(:,:,11) = -0.03;
constraint(:,:,12) = -0.0001;

constraint(:,:,13) = -0.03; %made it rlly small (prev 0003
constraint(:,:,14) = -0.0001;

constraint(:,:,15) = 0.003;
constraint(:,:,16) = 0.025;

constraint(:,:,17) = -0.01;
constraint(:,:,18) = 0.01;

constraint(:,:,19) = 0.003;
constraint(:,:,20) = 0.025;

constraint(:,:,21) = -0.01;
constraint(:,:,22) = 0.01;

constraint(:,:,23) = 0.98;
constraint(:,:,24) = 1.02;

constraint(:,:,25) = -0.002;
constraint(:,:,26) = 0.002;


[initialGuess, constraint]=global_guess_pf(pdf);
initialGuess = reshape(initialGuess, [imgPts, 13]); %%srinsrin originally 13
initialGuess = transpose(initialGuess);
constraint = reshape(constraint, [imgPts, 26]); %%srinsrin originally 13
initconstraint = single(transpose(constraint));
contype= int32(ones(1, 12) * 3); %%srin i.e. constrain both upper and lower bounds, see constrainttype.m

[parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
    model_id, initialGuess, initconstraint, contype, tolerance, max_n_iterations, [1,1,0,0,0,0,0,0,0,0,0,0,0].', EstimatorID.LSE, xValues);
disp('round 1')
[parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
    model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [0,0,1,0,0,0,0,0,0,0,0,0,0].', EstimatorID.LSE, xValues);
disp('round 2')
[parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
    model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [0,0,0,1,0,1,0,0,0,0,0,0,0].', EstimatorID.LSE, xValues);
disp('round 3')
 [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
     model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [1,1,0,0,0,0,0,0,0,0,0,0,0].', EstimatorID.LSE, xValues);
disp('round 4')
[parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
    model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [0,0,1,0,0,0,0,0,0,0,0,0,0].', EstimatorID.LSE, xValues);
disp('round 5')
 [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
     model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [0,0,0,0,1,0,1,0,0,0,0,0,0].', EstimatorID.LSE, xValues);
disp('round 6')
 [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
     model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [1,1,0,0,0,0,0,0,0,0,0,0,0].', EstimatorID.LSE, xValues);
disp('round 7')
[parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
    model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [1,1,1,0,0,0,0,0,0,0,0,0,0].', EstimatorID.LSE, xValues);
disp('round 8')

% tolerance=1e-21;
% max_n_iterations=10000;
 [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
     model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [1,1,1,1,1,1,1,1,0,1,0,0,0].', EstimatorID.LSE, xValues);
disp('round 9')

[parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
     model_id, parameters, initconstraint, contype, tolerance, max_n_iterations, [1,1,1,1,1,1,1,1,1,1,1,0,0].', EstimatorID.LSE, xValues);
disp('round 10')

% [parameters, states, chiSquares, n_iterations, ~] = gpufit_constrained(gpudata, [], ...
%      model_id, initialGuess, initconstraint, contype, tolerance, max_n_iterations, [1,1,1,0,0,0,0,0,0,0,0,0,0].', EstimatorID.LSE, xValues);

%prev allow gap dev to vary

% failed fits for pixel with extrem chiSquared or values outside of the
% measurement freq. range
% states(chiSquares > 5e-4) = 5;

fit = make_fit_struct(fit, initialPreGuess, initialGuess, parameters, states, ...
    chiSquares, n_iterations, nRes, sizeX, sizeY, kwargs.diamond);
fit.freq = freq;
fit.binSize = binSize;

msg = sprintf('%i: final GPU fitting complete in: %.1f s', nRes, toc(tStart));
logMsg('info',msg,1,0);

if numel(nonzeros(states)) > 0
    badPre = numel(nonzeros(states));
    msg = sprintf('%i: %i / %i pixels failed the final fit. See fit.states!', nRes, badPre, imgPts);
    logMsg('warn',msg,1,0);
end

if kwargs.checkPlot
    msg = sprintf('close figure to continue');
    logMsg('ATTENTION',msg,1,0);
    fig = gpu_fit_checkPlot(fit, binDataNorm, freq, binSize, kwargs.diamond);
    waitfor(fig)
end

end

%%
function data = slope_correction(data, freq, nPoints)
%[data] = slope_correction(data, freq, nPoints)
% calculates slope between 1st - last pixel and removes this from data
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
%[initialGuess] = get_initial_guess(gpudata, freq, diamond)
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
    case {'N15', 'doublet', 'gaussian', 'DAC', 'singlet', 'gauss4', 'gauss2', 'gaussTwin'}
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

function fit = make_fit_struct(fit, preGuess, initialGuess, parameters, states, chiSquares, n_iterations, nRes, sizeX, sizeY, diamond)
%[fit] = make_fit_struct(fit, preGuess, initialGuess, parameters, states, chiSquares, n_iterations, nRes, sizeX, sizeY, diamond)
msg = sprintf('%i: reshaping data into (%4i, %4i)', nRes, sizeY, sizeX);
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
fit.nRes = nRes;

if ~strcmp(preGuess, 'none')
    fit.pg = double(parameters_to_guess(preGuess, diamond));
else
    fit.pg = preGuess;
end

end
