function output = plot_peaks(params, peaks, freq, style)
    output = 0;
    switch style
        case "l"
            output = lorentzian(freq, peaks(1), params.A1, params.width);
            output = output + lorentzian(freq, peaks(2), params.A1, params.width);
            output = output + lorentzian(freq, peaks(3), params.A2, params.width);
            output = output + lorentzian(freq, peaks(4), params.A2, params.width);    
            output = output + lorentzian(freq, peaks(5), params.A3, params.width);
            output = output + lorentzian(freq, peaks(6), params.A3, params.width);    
            output = output + lorentzian(freq, peaks(7), params.A4, params.width);
            output = output + lorentzian(freq, peaks(8), params.A4, params.width);
        case "g"
            output = gaussian(freq, peaks(1), params.A1, params.width);
            output = output + gaussian(freq, peaks(2), params.A1, params.width);
            output = output + gaussian(freq, peaks(3), params.A2, params.width);
            output = output + gaussian(freq, peaks(4), params.A2, params.width);    
            output = output + gaussian(freq, peaks(5), params.A3, params.width);
            output = output + gaussian(freq, peaks(6), params.A3, params.width);    
            output = output + gaussian(freq, peaks(7), params.A4, params.width);
            output = output + gaussian(freq, peaks(8), params.A4, params.width);
        case "v"
            eta = 0.5;
            output = pseudovoight(freq, peaks(1), params.A1, params.width, eta);
            output = output + pseudovoight(freq, peaks(2), params.A1, params.width, eta);
            output = output + pseudovoight(freq, peaks(3), params.A2, params.width, eta);
            output = output + pseudovoight(freq, peaks(4), params.A2, params.width, eta);    
            output = output + pseudovoight(freq, peaks(5), params.A3, params.width, eta);
            output = output + pseudovoight(freq, peaks(6), params.A3, params.width, eta); 
            output = output + pseudovoight(freq, peaks(7), params.A4, params.width, eta);
            output = output + pseudovoight(freq, peaks(8), params.A4, params.width, eta);
    end
end

function output = lorentzian(x, x0, A, gamma)
    output = A/pi * ((gamma/2) ./ ((x - x0).^2 + (gamma/2).^2));
end

function output = gaussian(x, x0, A, sigma)
    output = A*exp(-(x-x0).^2 ./ (2 * sigma^2))/(sigma * sqrt(2*pi));
end

function output = pseudovoight(x, x0, A, sigma, eta)
    output = eta * gaussian(x, x0, A, sigma) + (1 - eta) * lorentzian(x, x0, A, sqrt(2 * log(2)) * sigma);
end

