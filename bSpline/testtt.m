% Knot vector and control points
knots = [0 0 0 1 2 3 4 4 4]; % Uniform B-spline
controlPoints = [0 1 2 3 2 1]; % Control points

% Create the B-spline
sp = spmak(knots, controlPoints);

% Evaluate the B-spline at specific points
x = linspace(0, 4, 100); % Evaluation points
y = fnval(sp, x); % Evaluate the spline

% Plot
plot(x, y, 'b', controlPoints, 'ro-');
legend('B-spline curve', 'Control points');
xlabel('x');
ylabel('y');
title('B-spline');
grid on;
