function find_line_cuts(data, val, dir)
    if strcmp(dir, "Vertical")
        data=squeeze(data(val, :, :));
    elseif strcmp(dir, "Horizontal")
        data=squeeze(data(:, val, :));
    else
        data=squeeze(data(:, val, :));
    end
    % Normalize data so linecuts can actually be taken
    for i=1:size(data,1)
        data(i,:)=data(i,:)-mean(data(i,:));
        data(i,:)=data(i,:)/std(data(i,:));
    end
    for xsmooth=1:size(data, 1)
        data(xsmooth, :) = smooth(data(xsmooth, :), 41, 'sgolay', 7);
    end
    % data = data;

    % Take data and convert it to a black and white image
    figure;
    subplot(1, 4, 1);    
    imagesc(data);
    ylim([1 size(data,1)]);
    ylabel('Pixel');
    ax = gca;
    ax.YDir = 'reverse';

    subplot(1, 4, 2);
    T = imbinarize(data, adaptthresh(data, 0.5, "ForegroundPolarity","bright"));
    imagesc(T);
    ylim([1 size(data,1)]);
    ylabel('Pixel');
    ax = gca;
    ax.YDir = 'reverse';

    subplot(1, 4, 3);
    T = medfilt2(T);
    % T = bwareaopen(T, 1000, 8);
    % T = imfill(T, 'holes');
    % T = imcomplement(T);
    % T = bwareaopen(T, 1000, 8);
    % T = imfill(T, 'holes');
    % T = imcomplement(T);
    imagesc(T);
    ylim([1 size(data,1)]);
    ylabel('Pixel');
    ax = gca;
    ax.YDir = 'reverse';

    subplot(1, 4, 4);
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

