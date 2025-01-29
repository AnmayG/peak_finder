function [guess, constraint] = global_guess_pf(pdf, kwargs)
%[guess] = global_guess(data, freq; 'forceGuess', 'checkPlot', 'smoothDegree', 'minPeakDistance')
% Returns a global guess for the given dataset
%
% Parameters
% ----------
% 
% Returns
% -------
%     array with 6 x X x Y elements, where X and Y are the indices of the
%     pixels

arguments
    pdf
    % keyword arguments
    kwargs.forceGuess (1,1) {mustBeBoolean(kwargs.forceGuess)} = false
    kwargs.checkPlot (1,1) {mustBeBoolean(kwargs.checkPlot)} = false
    kwargs.smoothDegree  (1,1) {mustBeNumeric, mustBePositive} = 2
    kwargs.minPeakDistance (1,1) {mustBeNumeric} = 0
end

%%
msg = sprintf('generating initial guess from global resonance parameters');
logMsg('info',msg,1,0);

sizeX = size(pdf,2); % binned image x-dimensions
sizeY = size(pdf,1); % binned image y-dimensions
% meanData = squeeze(mean(mean(data,1),2));
% m1 = mean(data,1);
% m2 = nanmean(squeeze(m1),1);
meanData = squeeze(mean(pdf,[1,2],'omitnan'));

%%srin I'm commenting out the below
%% Resonance 1:
% [pkVal, pkLoc, fitFlg] = guess_peaks(meanData, meanData, freq, ...
%                                'smoothDegree', kwargs.smoothDegree, ...
%                                'forceGuess', kwargs.forceGuess,...
%                                'checkPlot', kwargs.checkPlot);
% 
% Rguess = [(pkLoc(1)+pkLoc(2)+pkLoc(3))/3  0.0005  ( mean(meanData(1:10)) +pkVal-1)' mean(meanData(1:10))-1 ]; %resonance [GHz], Width [GHz], (contrast1; contrast2; contrast3)', baseline
% Rguess = Rguess';
% 
% % write guess
% guess = zeros(sizeY, sizeX, 6);
% guess(:,:,1) = Rguess(1);
% guess(:,:,2) = Rguess(2);
% guess(:,:,3) = Rguess(3);
% guess(:,:,4) = Rguess(4);
% guess(:,:,5) = Rguess(5);
% guess(:,:,6) = Rguess(6);
% 
% guess = single(guess);


%%srin first by hand I'm just going to initialize all 13 fit parameters as
%%the same everywhere

% * parameters: An input vector of model parameters.
% *             p[0]: D111
% *             p[1]: E111
% *             p[2]: std 111 avg
% *             p[3]: std 111 dev
% *             p[4]: contrast 111 avg
% *             p[5]: contrast 111 dev
% *             p[6]: Dnon111 dev (respect to D111)
% *             p[7]: Enon111 dev (respect to E111/3 i.e. z field)
% *             p[8]: std non111 avg
% *             p[9]: std non111 dev
% *             p[10]: contrast non111 avg
% *             p[11]: contrast non111 dev
% *             p[12]: baseline


%%in matlab



% *             p[1]: D111, average of 111 peaks location
% *             p[2]: E111, half the difference of 111 peaks location
% *             p[3]: std 111 avg, average standard dev of the 111 peaks
% *             p[4]: std 111 dev, difference in the standard dev of the 111 peaks
% *             p[5]: contrast 111 avg, average height of the 111 peaks 
% *             p[6]: contrast 111 dev, difference in the heights of the 111 peaks
% *             p[7]: Dnon111 dev, difference of the center of the non111s (with respect to D111)
% *             p[8]: Enon111 dev, difference in the E of the non 111s (respect to E111/3 i.e. z field)
% *             p[9]: std non111 avg, average standard dev of the non111
% peaks DON'T MAKE THIS ZERO OR IT WILL DIVERGE
% *             p[10]: std non111 dev, difference in the standard dev of the non 111 peaks
% *             p[11]: contrast non111 avg, average height of the non 111 peaks 
% *             p[12]: contrast non111 dev, difference in the heights of the 111 peaks
% *             p[13]: baseline

%%constraints should be 
%1: within some pressure std (peakfind)
% 2: within some field std (peakfind)
% 3: some width std
% 4: from -'3' to +'3', starting at 0
% 5: some height std
% 6: 

%2409 100G 16.5k
guess = zeros(sizeY, sizeX, 13);

%[1,1,1,0,1,0,0,0,0,0,0,0,1] is the variation array
% guess(:,:,1) = 4.238;
% guess(:,:,2) = 0.274;
% guess(:,:,3) = 0.062;
% guess(:,:,4) = -0.034;
% guess(:,:,5) = -0.013; %-0.01;
% guess(:,:,6) = 0.0102;
% guess(:,:,7) = 0;
% guess(:,:,8) = 0;
% guess(:,:,9) = 1;
% guess(:,:,10) = 0;
% guess(:,:,11) = 0; %-0.01;
% guess(:,:,12) = 0;
% guess(:,:,13) = 1;


% %2409 100g 285k
% 
% guess(:,:,1) = 4.137;
% guess(:,:,2) = 0.284;
% guess(:,:,3) = 0.062;
% guess(:,:,4) = -0.034;
% guess(:,:,5) = -0.013; %-0.01;
% guess(:,:,6) = 0.0102;
% guess(:,:,7) = 0;
% guess(:,:,8) = 0;
% guess(:,:,9) = 1;
% guess(:,:,10) = 0;
% guess(:,:,11) = 0; %-0.01;
% guess(:,:,12) = 0;
% guess(:,:,13) = 1;


% %srin 2403 40G/50G
% guess = zeros(sizeY, sizeX, 13);
% 
% guess(:,:,1) = 2.92;
% guess(:,:,2) = 0.0028*50;
% guess(:,:,3) = 0.010;
% guess(:,:,4) = 0;
% guess(:,:,5) = -0.01; %-0.01;
% guess(:,:,6) = 0;
% guess(:,:,7) = 0;
% guess(:,:,8) = 0;
% guess(:,:,9) = 0.010;
% guess(:,:,10) = 0;
% guess(:,:,11) = -0.01; %-0.01;
% guess(:,:,12) = 0;
% guess(:,:,13) = 1; %%srin seems like it normalizes in prepare_raw_data

% %srin nickelate from z/e
% guess = zeros(sizeY, sizeX, 13);
% 
% guess(:,:,1) = 3.08;
% guess(:,:,2) = 0.0028*150;
% guess(:,:,3) = 0.008;
% guess(:,:,4) = 0.002;
% guess(:,:,5) = -0.007; %-0.01;
% guess(:,:,6) = -0.001;
% guess(:,:,7) = 0;
% guess(:,:,8) = 0;
% guess(:,:,9) = 1;
% guess(:,:,10) = 0;
% guess(:,:,11) = 0; %-0.01;
% guess(:,:,12) = 0;
% guess(:,:,13) = 1; %%srin seems like it normalizes in prepare_raw_data

%quadratic
% guess = zeros(sizeY, sizeX, 3);
% 
% guess(:,:,1) = 0.5;
% guess(:,:,2) = 1;
% guess(:,:,3) = 0;



%%2416 16gpa

% 
% guess = zeros(sizeY, sizeX, 13);
% 
% guess(:,:,1) = 3.05;
% guess(:,:,2) = 0.0028*206.6;
% guess(:,:,3) = 0.020; %0.01
% guess(:,:,4) = 0;
% guess(:,:,5) = -0.01; %-0.01;
% guess(:,:,6) = 0;
% guess(:,:,7) = 0;
% guess(:,:,8) = 0;
% guess(:,:,9) = 0.010;
% guess(:,:,10) = 0;
% guess(:,:,11) = -0.01; %-0.01;
% guess(:,:,12) = 0;
% guess(:,:,13) = 1; %%srin seems like it normalizes in prepare_raw_data
% 
% 

% guess = zeros(sizeY, sizeX, 13);
% 
% guess(:,:,1) = 3.1;
% guess(:,:,2) = 0.0028*5;
% guess(:,:,3) = 0.020; %0.01
% guess(:,:,4) = 0;
% guess(:,:,5) = -0.01; %-0.01;
% guess(:,:,6) = 0;
% guess(:,:,7) = 0;
% guess(:,:,8) = 0;
% guess(:,:,9) = 0.010;
% guess(:,:,10) = 0;
% guess(:,:,11) = -0.01; %-0.01;
% guess(:,:,12) = 0;
% guess(:,:,13) = 1; %%srin seems like it normalizes in prepare_raw_data
% 
% guess = zeros(sizeY, sizeX, 12);
% 
% guess(:,:,1) = 3.1;
% guess(:,:,2) = 0.0028*206.6;
% guess(:,:,3) = 0.5*25e-3; %0.01
% guess(:,:,4) = -0.01;
% guess(:,:,5) = 0; %-0.01;
% guess(:,:,6) = -0.01;
% guess(:,:,7) = 0;
% guess(:,:,8) = 0.01;
% guess(:,:,9) = 0;
% guess(:,:,10) = 0.01;
% guess(:,:,11) = 0; %-0.01;
% guess(:,:,12) = 1;

%put in weights that's just value
% guess = zeros(sizeY, sizeX, 13);
% 
% guess(:,:,1) = 3.1;
% guess(:,:,2) = 0.0028*206.6;
% guess(:,:,3) = -0.2*25e-3; %0.01
% guess(:,:,4) = -0.01;
% guess(:,:,5) = -0.005; %-0.01;
% guess(:,:,6) = -0.01;
% guess(:,:,7) = -0.005;
% guess(:,:,8) = 0.008;
% guess(:,:,9) = 0;
% guess(:,:,10) = 0.008;
% guess(:,:,11) = 0; %-0.01;
% guess(:,:,12) = 1;
% guess(:,:,13) = 0;


% guess = zeros(sizeY, sizeX, 13);
% 
% guess(:,:,1) = 3.05;
% guess(:,:,2) = 0.203; %0.0028*206.6;
% guess(:,:,3) = -0.01;%-0.05*25e-3; %0.01
% 
% guess(:,:,4) = -0.01;
% guess(:,:,5) = -0.005; %-0.01;
% guess(:,:,6) = -0.01;
% guess(:,:,7) = -0.005;
% 
% guess(:,:,8) = 0.008;
% guess(:,:,9) = 0;
% guess(:,:,10) = 0.008;
% guess(:,:,11) = 0; %-0.01;
% guess(:,:,12) = 1;
% guess(:,:,13) = 0;

%pkf=all_param_data.parameters_dataframe;
% apd = evalin('base', 'all_param_data'); %things in the workspace need to be called in explicitly
% pkf=apd.parameters_dataframe;
guess = zeros(sizeY, sizeX, 13);
guess(:,:,1)=pdf(:,:,1)/1e9;
guess(:,:,2)=pdf(:,:,2)/1e9;
guess(:,:,3)=0.001; %prev -0.01

guess(:,:,4)=1/3*(pdf(:,:,5)/1e9-pdf(:,:,6)/1e9);
guess(:,:,5)=2/3*(pdf(:,:,5)/1e9-pdf(:,:,6)/1e9);
guess(:,:,6)=1/3*(pdf(:,:,5)/1e9+pdf(:,:,6)/1e9);
guess(:,:,7)=2/3*(pdf(:,:,5)/1e9+pdf(:,:,6)/1e9);

guess(:,:,8)=pdf(:,:,3)/1e9; %bc peakfinder already converts FWHM to std
guess(:,:,9)=pdf(:,:,4)/1e9;
guess(:,:,10)=pdf(:,:,3)/1e9; 
guess(:,:,11)=pdf(:,:,4)/1e9;

guess(:,:,12) = 1;
guess(:,:,13) = 0;


constraint = zeros(sizeY, sizeX, 26);
constraint(:,:,1) = pdf(:,:,1)/1e9-0.1;
constraint(:,:,2) = pdf(:,:,1)/1e9+0.1;

constraint(:,:,3) = pdf(:,:,2)/1e9-0.1;
constraint(:,:,4) = pdf(:,:,2)/1e9+0.1;

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


%the below line is essential
guess = single(guess);

