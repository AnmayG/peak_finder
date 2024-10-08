function [transForm, refFrame] = get_image_tform_complex(fixedData, movingData, kwargs)
%[transForm, refFrame] = get_image_tform_complex(fixedData, movingData; 'checkPlot', 'transformationType', 'binning')
% Function lets you pick several points on a reference image and the target
% image. It calculates a transformation
% 
% Parameters
% ----------
%     fixedData: double QDM/LED data
%     movingData: double QDM/LED data
%         data to be matched to the refernce data
%     binning: int [1]
%       needed for correct transformation
%     checkPlot: bool [false]
%         Adds a plot to check alignment if true
%
% Note
% ----
%     Both fixedData and movingData have to be data, not jpg files. You can
%     load a jpg with imread().
arguments
    fixedData double
    movingData double
    kwargs.checkPlot {mustBeBoolean(kwargs.checkPlot)} = false 
    kwargs.transformationType {mustBeMember(kwargs.transformationType , ['similarity','projective','polynomial'])} = 'similarity' 
    kwargs.binning {mustBeInteger(kwargs.binning)} = 1
end 

if max(fixedData, [], 'all') > 1
    fixedData = double(fixedData);
    fixedData = fixedData / max(fixedData, [], 'all');
end

if max(movingData, [], 'all') > 1
    movingData = double(movingData);
    movingData = movingData / max(movingData, [], 'all');
end

[mp, fp] = cpselect(movingData, fixedData, 'Wait', true);

if isequal(kwargs.transformationType, 'polynomial')
    transForm = fitgeotrans(mp, fp, kwargs.transformationType, 2);
else
    transForm = fitgeotrans(mp, fp, kwargs.transformationType);
end

% create refence frame for fixed image
refFrame = imref2d(size(fixedData));

if isequal(kwargs.checkPlot, true)
    checkFigure = figure('Name', 'Align images');
    
    transformedData = imwarp(movingData, transForm, 'OutputView', refFrame);

    subplot(2, 1, 1)
    imshowpair(fixedData, movingData, 'Scaling', 'joint')
    title('uncorrected')
    axis xy; axis equal, axis tight

    subplot(2, 1, 2)
    imshowpair(fixedData, transformedData, 'Scaling', 'joint')
    title('corrected')
    axis xy; axis equal, axis tight

end