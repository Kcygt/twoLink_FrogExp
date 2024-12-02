% Sigmoid time function and its derivative in MATLAB

% Parameters
T = 10; % Total time duration (seconds)
x0 = T / 2; % Midpoint of the transition
k = 5; % Steepness of the transition

% Time vector (from 0 to T)
t = linspace(0, T, 1000);

% Sigmoid function
sigmoid = 1 ./ (1 + exp(-k * (t - x0)));

% Derivative of the sigmoid function
sigmoid_derivative = (k * exp(-k * (t - x0))) ./ (1 + exp(-k * (t - x0))).^2;

% Plot the sigmoid function and its derivative
figure;

subplot(2,1,1);
plot(t, sigmoid, 'LineWidth', 2);
title('Sigmoid Time Function');
xlabel('Time (seconds)');
ylabel('Sigmoid Value');
grid on;
xlim([0 T]);
ylim([0 1]);

subplot(2,1,2);
plot(t, sigmoid_derivative, 'LineWidth', 2);
title('Derivative of Sigmoid Time Function');
xlabel('Time (seconds)');
ylabel('Derivative');
grid on;
xlim([0 T]);

