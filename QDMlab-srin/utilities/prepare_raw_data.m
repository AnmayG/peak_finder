function [binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes, header, kwargs)
%[binDataNorm, freq] = prepare_raw_data(expData, binSize, nRes; 'header', 'gpuData', 'crop', 'normalize', 'fcrop')
% prepares the raw data for GPU fitting
% 1. reshapes the data into from (x*y) -> (y,x) array
% 2. Bins data: (imresize)
% 3. Normalizes
%%srin this thing needs to be fed 'gWide' directly
% Parameters
% ----------
%     required
%     ========
%     expData: struct
%         Data of load(run0000n.mat)
%     binSize: int
%         binning size (can be 1)
%     header: ('none')
%       in RARE cases a header file is needed, if the frequency data is not
%       stored in the data file 
%     nRes: int
%         number of resonance. Low frequencies = 1, High frequencies = 2
%     crop: (int, int, int, int) ['none']
%         cell with {x0, y0, dx, dy} where x0 and y0 is the lower left
%         corner. See also: pick_box
%         if crop = 'none' -> no cropping
%     normalize: bool [true]
%         Chose if the data should be normalized

arguments
    expData
    binSize = 1 %%srin
    nRes = 1 %%srin
    header = 'none'
    kwargs.gpuData = false
    kwargs.crop = 'none'
    kwargs.normalize = true
    kwargs.fcrop  = false
end
%srin freq is the list of frequencies, because they do 2 sweeps they have to
%make two lists. instead, just pull directly from gWide

%%dataStack = expData.(sprintf('imgStack%i',nRes));

%%fRanges = get_franges(expData, header);   %everything is in GHz

%%freq = fRanges{nRes};
freq=expData.SweepParam/1E9; %%srin let's use GHz since it makes values like E and deltaD nicer

%% data preparation
% X/Y of unbinned data
% Note: X = COL; Y = ROW -> (y,x) for (row, col) matlab convention
%%spanXTrans = 1:expData.imgNumCols;
%%spanYTrans = 1:expData.imgNumRows;

spanXTrans = 1:size(expData.signal,1);
spanYTrans = 1:size(expData.signal,2);

%%srin this thing is trying to stack the low and high mw frequencies, which
%%we don't need to do. even so, retain datastack as a variable for the
%%cropping
% check for 101 frequencies. File includes imgStack3
% % if isfield(expData, 'imgStack3')      
% %     % combine 1&2 or 3&4
% %     if nRes == 1
% %         a = 1; b = 2;
% %     else
% %         a = 3; b = 4;
% %     end
% %     dataStacka = expData.(sprintf('imgStack%i',a)); 
% %     dataStackb = expData.(sprintf('imgStack%i',b));
% %     dataStack = [dataStacka; dataStackb];
% % end

if kwargs.fcrop
    dataStack = dataStack(kwargs.fcrop(1):kwargs.fcrop(2), :);
    freq = freq(kwargs.fcrop(1):kwargs.fcrop(2));
end
% reshape and transpose each image
% NOTE: I guess the data is written pixel by pixel, starting left top
% then row by row -> freq x col * rows (i.e. [51, 2304000])

%%srin their data is natively in a line, but ours is in a grid already, so
%%we don't need to reshape here

%%data = QDMreshape(dataStack, expData.imgNumRows, expData.imgNumCols);
data=expData.signal;

%% Data cropping
if ~strcmp(kwargs.crop, 'none') 
    % cropping
    x0 = kwargs.crop(1);
    y0 = kwargs.crop(2);
    x1 = kwargs.crop(3)+x0;
    y1 = kwargs.crop(4)+y0;
    
    msg = sprintf('<>   %i: cropping data >> crop = X(%i, %i) Y(%i, %i)', ...
        nRes, x0,x1,y0,y1);
    logMsg('debug',msg,1,0);
    data = data(y0:y1,x0:x1,:);
end
   
%% binning
msg = sprintf('<>   %i: binning data >> binSize = %i', nRes, binSize);
logMsg('debug',msg,1,0);
binData = imresize(data, 1/binSize, 'method', 'box');

% Correct for severely non-unity baseline by dividing pixelwise by
% average of all frequency points
if kwargs.normalize
    binDataNorm = zeros(size(binData));
    NormalizationFactor = mean(binData,3);    % compute average
    for y = 1:length(freq)
        binDataNorm(:,:,y) = binData(:,:,y) ./ NormalizationFactor;
    end
else
    binDataNorm = binData;
end

%% return gpudata if kwargs.gpuData == true
if isequal(kwargs.gpuData, true)
    binDataNorm = gpu_data_reshape(binDataNorm);
end

