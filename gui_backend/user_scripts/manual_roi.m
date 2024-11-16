function function_outputs = manual_roi(function_inputs)
    %   function_inputs = struct array with fields:
    %       signal: [n × n × freq double]
    %   SweepParam: freq
    %   parameters: [n × n × 17 double]
    %    peak_info: [1×1 struct] with fields locs, vals, widths, proms
    %   reflection: [n × n double]
    %     settings: [1×1 struct] with each parameter from inputs
    % new_settings: [1×1 struct] with parameters at time of run
    function_outputs = function_inputs;
    % Modify function_outputs fields below
    prev_params = function_inputs.parameters;
    new_mask = paintMask(prev_params(:, :, 1));
    for i=length(function_inputs)
        function_outputs(i).settings(i).roi_mask = new_mask;
    end
    function_outputs.new_settings.roi_mask = new_mask;
end

function mask = paintMask(data)
    figureHandle = figure('Name', 'Paint Mask -- Left Click to Paint, Right Click to Erase');
    data(data == 0) = NaN;
    matrix = data;
    imagesc(data, 'AlphaData', ~isnan(data));
    dataMin = min(matrix(:), [], 'omitnan');
    dataMax = max(matrix(:), [], 'omitnan');
    clim([dataMin dataMax]);
    axis equal;
    axis tight;
    colorbar;
    hold on;
    mask = zeros(size(matrix));
    
    brushSize = 5;
    halfBrush = floor(brushSize / 2);
    
    set(figureHandle, 'WindowButtonDownFcn', @mouseDown);
    set(figureHandle, 'WindowButtonUpFcn', @mouseUp);
    set(figureHandle, 'CloseRequestFcn', @closeFigure);

    isPainting = false;
    isErasing = false;
    maskOverlay = imagesc(mask);
    set(maskOverlay, 'AlphaData', 0); % Initially invisible

    function mouseDown(~, ~)
        % Determine if painting or erasing based on mouse button
        clickType = get(figureHandle, 'SelectionType');
        if strcmp(clickType, 'normal')  % Left-click for painting
            isPainting = true;
            set(figureHandle, 'WindowButtonMotionFcn', @paint);
            paint();
        elseif strcmp(clickType, 'alt')  % Right-click for erasing
            isErasing = true;
            set(figureHandle, 'WindowButtonMotionFcn', @erase);
            erase();
        end
    end

    function mouseUp(~, ~)
        isPainting = false;
        isErasing = false;
        set(figureHandle, 'WindowButtonMotionFcn', '');
    end

    function paint(~, ~)
        if isPainting
            currentPoint = get(gca, 'CurrentPoint');
            x = round(currentPoint(1, 1));
            y = round(currentPoint(1, 2));
            
            if x > 0 && y > 0 && x <= size(matrix, 2) && y <= size(matrix, 1)
                xRange = max(1, x-halfBrush):min(size(matrix, 2), x+halfBrush);
                yRange = max(1, y-halfBrush):min(size(matrix, 1), y+halfBrush);            
                mask(yRange, xRange) = 1;
                set(maskOverlay, 'CData', mask);
                set(maskOverlay, 'AlphaData', mask * 0.5); % Semi-transparent mask
            end
        end
    end

    function erase(~, ~)
        if isErasing
            currentPoint = get(gca, 'CurrentPoint');
            x = round(currentPoint(1, 1));
            y = round(currentPoint(1, 2));

            if x > 0 && y > 0 && x <= size(matrix, 2) && y <= size(matrix, 1)
                xRange = max(1, x-halfBrush):min(size(matrix, 2), x+halfBrush);
                yRange = max(1, y-halfBrush):min(size(matrix, 1), y+halfBrush);
                
                mask(yRange, xRange) = 0;
                set(maskOverlay, 'CData', mask);
                set(maskOverlay, 'AlphaData', mask * 0.5); % Semi-transparent mask
            end
        end
    end

    function closeFigure(~, ~)
        delete(figureHandle);
    end
    waitfor(figureHandle);
end




