function [transform, refFrame, success] = align_images(fixed_data, moving_data, display)
    % Template code from
    % https://www.mathworks.com/help/vision/ref/matchfeatures.html
    success = true;
    fixed_data = rescale(fixed_data);
    moving_data = rescale(moving_data);
    transform = affine2d([[1 0 0]; [0 1 0]; [0 0 1]]);
    refFrame = imref2d(size(fixed_data));

    if isequal(fixed_data, moving_data)
        % Transformation not needed since duplicate image
        return
    end

    % KAZE features did a slight rotation which was good
    % ORB features also moved it down a bit
    % 
    ptsOriginal  = detectKAZEFeatures(fixed_data);
    ptsDistorted = detectKAZEFeatures(moving_data);
    
    [featuresOriginal,validPtsOriginal] = extractFeatures(fixed_data,ptsOriginal);
    [featuresDistorted,validPtsDistorted] = extractFeatures(moving_data,ptsDistorted);

    index_pairs = matchFeatures(featuresOriginal, featuresDistorted);

    matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
    matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));

    try
        [transform, inlierIdx] = ...    
            estgeotform2d(matchedPtsDistorted,matchedPtsOriginal,...    
                'affine');
    catch ME
        report = getReport(ME);
        disp(report)
        display_images(fixed_data, moving_data, transform, refFrame, ...
            ptsOriginal, ptsDistorted)
        success = false;
        return;
    end

    refFrame = imref2d(size(fixed_data));
    inlierPtsDistorted = matchedPtsDistorted(inlierIdx,:);
    inlierPtsOriginal  = matchedPtsOriginal(inlierIdx,:);


    if display
        display_images(fixed_data, moving_data, transform, refFrame, ...
            ptsOriginal, ptsDistorted, inlierPtsOriginal, inlierPtsDistorted)
    end
end

function display_images(fixed_data, moving_data, transform, refFrame, ...
    ptsOriginal, ptsDistorted, inlierPtsOriginal, inlierPtsDistorted)
    plot_inliers = true;
    if nargin < 7
        plot_inliers = false;
    end
    % It's time to do some sketchy boolean casting
    transformedData = imwarp(moving_data,transform, 'OutputView', refFrame);    
    figure;    
    
    subplot(1 + plot_inliers, 3, 1)
    imshowpair(fixed_data, moving_data, 'Scaling', 'joint')    
    title('uncorrected')    
    axis xy; 
    axis equal; 
    axis tight;    
    a = gca;    
    
    if plot_inliers
        subplot(2, 3, 2)
        showMatchedFeatures(fixed_data,moving_data, inlierPtsOriginal,inlierPtsDistorted);    
        title('Matched features');    
        axis xy; 
        axis equal; 
        axis tight;    
        b = gca;    
        
        subplot(2, 3, 3)    
        imshowpair(fixed_data, transformedData, 'Scaling', 'joint')    
        title('corrected')    
        axis xy; 
        axis equal; 
        axis tight;    
        c = gca;    
        linkaxes([a b c])
        
        disp(invert(transform))
        subplot(2, 3, 6)
        bar3(inv(transform.A) * transform.A' - eye(3))
        title("Affine Transformation Matrix - Goodness: " + det( ...
            transform.A))
    end

    subplot(1 + plot_inliers, 3, 2 + plot_inliers * 2)
    imshow(fixed_data); 
    hold on;
    plot(ptsOriginal.selectStrongest(10));
    title("Original Features (Top 10 Strongest)")
    hold off;

    subplot(1 + plot_inliers, 3, 3 + plot_inliers * 2)
    imshow(moving_data); 
    hold on;
    plot(ptsDistorted.selectStrongest(10));
    title("Distorted Features (Top 10 Strongest)")
    hold off;
end
