function [outputArg1,outputArg2] = find_line_cuts(data, val, dir)
    if strcmp(dir, "Vertical")
        data=squeeze(data(val, :, :));
    elseif strcmp(dir, "Horizontal")
        data=squeeze(data(:, val, :));
    else
        data=squeeze(data(:, val, :));
    end
    % Normalize data
    for i=1:size(data,1)
        data(i,:)=data(i,:)-mean(data(i,:));
        data(i,:)=data(i,:)/std(data(i,:));
    end
    % data(data > 0) = 0;
    % data = normalize(data);
    data = data;

    % Take data and convert it to a black and white image
    figure;
    subplot(1, 3, 2);
    imagesc(imbinarize(data, adaptthresh(data, 0.5, "ForegroundPolarity","dark")));
    ylim([1 size(data,1)]);
    ylabel('Pixel');
    ax = gca;
    ax.YDir = 'reverse';

    subplot(1, 3, 1);    
    imagesc(data);
    ylim([1 size(data,1)]);
    ylabel('Pixel');
    ax = gca;
    ax.YDir = 'reverse';
    colorbar;

    subplot(1, 3, 3);
    colormap gray;
    disp(graythresh(data))
    T = edge(data, 'log');
    T = bwareafilt(T, [10, Inf], 8);
    hold on;
    axis tight;
    data(T) = NaN;
    imagesc(data);
    % imagesc(T);
    hold off;
    ylim([1 size(T,1)]);
    ylabel('Pixel');
    ax = gca;
    ax.YDir = 'reverse';
end

