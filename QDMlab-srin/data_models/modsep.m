function modsep(p, x)
%plots gausstwinD but with the four peaks separated, call like 
%modsep(squeeze(f8.parameters(:,yy,xx)),freq)
    % Ensure p is the correct size
    if length(p) ~= 13
        error('Input parameter p must be a vector of length 13.');
    end


    % Calculate Gaussian components
    eps = 1e-6; % Small epsilon to avoid division by zero or negative variances
    p(8) = max(p(8), eps);
    p(9) = max(p(9), eps);
    p(10) = max(p(10), eps);
    p(11) = max(p(11), eps);

    argx1 = (x - (p(1) - p(2) - p(3) + p(13))).^2 ./ (2 * (p(8) - p(9)).^2);
    argx2 = (x - (p(1) - p(2) + p(3) - p(13))).^2 ./ (2 * (p(10) - p(11)).^2);
    argx3 = (x - (p(1) + p(2) - p(3) - p(13))).^2 ./ (2 * (p(8) + p(9)).^2);
    argx4 = (x - (p(1) + p(2) + p(3) + p(13))).^2 ./ (2 * (p(10) + p(11)).^2);

    ex1 = exp(-argx1);
    ex2 = exp(-argx2);
    ex3 = exp(-argx3);
    ex4 = exp(-argx4);

    % Weighted components
    g1 = p(4) * ex1;
    g2 = p(5) * ex2;
    g3 = p(6) * ex3;
    g4 = p(7) * ex4;

    % Sum of all components
    model = g1 + g2 + g3 + g4 + p(12);

    % Plot the Gaussians and their sum
    figure;
    hold on;
    plot(x, g1+ p(12), 'LineWidth', 2, 'DisplayName', 'Gaussian 1');
    plot(x, g2+ p(12), 'LineWidth', 2, 'DisplayName', 'Gaussian 2');
    plot(x, g3+ p(12), 'LineWidth', 2, 'DisplayName', 'Gaussian 3');
    plot(x, g4+ p(12), 'LineWidth', 2, 'DisplayName', 'Gaussian 4');
    plot(x, model, '--', 'LineWidth', 2, 'DisplayName', 'Sum of Gaussians');

    % Add labels, legend, and grid
    title('Gaussian Components and their Sum');
    xlabel('x');
    ylabel('Amplitude');
    legend;
    grid on;
    hold off;
end
