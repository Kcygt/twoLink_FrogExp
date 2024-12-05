% Given points
x = [0 1 2 3];
y = [0 1 1 0];

% Define the knot vector for B-spline
% A clamped B-spline requires adding repeated knots at the start and end
knots = [0 0 0 1 2 3 3 3]; % Cubic B-spline (degree 3)

% Create a B-spline using spapi (spline approximation)
degree = 3; % Degree of the spline
spline_function = spapi(knots, x, y);

% Plot the original points
figure;
plot(x, y, 'o', 'MarkerSize', 10, 'DisplayName', 'Original Points');
hold on;

% Evaluate and plot the B-spline
xx = linspace(min(x), max(x), 100);
yy = fnval(spline_function, xx);
plot(xx, yy, '-r', 'LineWidth', 1.5, 'DisplayName', 'B-spline Curve');

% Add labels and legend
xlabel('x');
ylabel('y');
title('B-Spline Curve');
legend show;
grid on;
hold off;
