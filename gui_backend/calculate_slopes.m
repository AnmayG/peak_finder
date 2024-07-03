function [slope, r_squared, intercept] = calculate_slopes(y_vals, x_coords, zfs)
        X = [ones(length(x_coords), 1) x_coords];
        % Effect = @(b,x)(2.8025*sqrt(zfs^2 + abs(b^2 * x.^2)));
        opts = optimset('Display','off');
        fitResult = lsqlin(X, y_vals, [], [], [], [], [0, 0], [inf, inf], [], opts);
        % disp(fitResult)
        % coefficients = X \ y_vals;
        slope = fitResult(2); % coefficients(2);
        intercept = fitResult(1); % coefficients(1);
        
        % Calculate r^2
        y_calc = X * fitResult;
        r_squared = 1 - sum((y_vals - y_calc).^2)/sum((y_vals - mean(y_vals)).^2);
end

