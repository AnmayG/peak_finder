function [baseline_value] = find_baseline(y, x, data, freq, threshold, display)
    if nargin < 5
        threshold = 0.4;
        display = false;
    elseif nargin < 6
        display = false;
    end
    cropped_signal = squeeze(data(x, y, :));
    peak_threshold = prctile(cropped_signal, threshold * 100);
    cropped_signal(cropped_signal < peak_threshold) = peak_threshold;
    baseline_value = mean(cropped_signal(cropped_signal ~= peak_threshold));
    
    if display
        figure
        hold on
        plot(gWide.SweepParam, cropped_signal, freq, squeeze(gWide.signal(x, y, :)))
        yline(peak_threshold, 'Color', 'g')
        yline(baseline_value, 'Color', 'b')
    end
end
