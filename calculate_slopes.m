function [slope, r_squared, intercept] = calculate_slopes(y_vals, x_coords)
%CALCULATE_SLOPES Calculates slopes, intercepts, and r^2 values for given
% y_vals and x_coords
%   If implementing for all_data, make sure to iterate through every (x, y)
        X = [ones(length(x_coords), 1) x_coords];
        coefficients = X \ y_vals;
        slope = coefficients(2);
        intercept = coefficients(1);
        
        % Calculate r^2
        y_calc = X * coefficients;
        r_squared = 1 - sum((y_vals - y_calc).^2)/sum((y_vals - mean(y_vals)).^2);
end

