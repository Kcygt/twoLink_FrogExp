% Parameters
K = .6; % Scalar value for K

% Define range for e
e = linspace(-30, 30, 100); % 100 points between -2 and 2

% Calculate the function
f = K * e .* exp(-e .* K .* e);
f = 1 ./ (1 + exp(-e));
% Plot the result
figure;
plot(e, f, 'LineWidth', 2);
xlabel('e');
ylabel('f(e)');
title('Plot of f(e) = K \cdot e \cdot exp(-e'' \cdot K \cdot e)');
grid on;
