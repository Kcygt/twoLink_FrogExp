% B-Spline Curve Generator
% Define 4 control points
control_points = [0, 0; 1, 2; 3, 3; 4, 0]; % Modify these points as needed

% Degree of the B-spline
degree = 3; % Cubic B-spline

% Number of points for the curve
num_points = 100;

% Generate knot vector
n = size(control_points, 1); % Number of control points
knot_vector = [zeros(1, degree+1), linspace(0, 1, n-degree), ones(1, degree+1)];

% Parameter values for the curve
t = linspace(knot_vector(degree+1), knot_vector(end-degree), num_points);

% Evaluate the B-spline curve
curve = zeros(num_points, 2);
for i = 1:num_points
    curve(i, :) = de_boor(control_points, knot_vector, degree, t(i));
end

% Plot the control points and the curve
figure;
hold on;
plot(control_points(:, 1), control_points(:, 2), 'ro-', 'LineWidth', 2, 'DisplayName', 'Control Points');
plot(curve(:, 1), curve(:, 2), 'b-', 'LineWidth', 2, 'DisplayName', 'B-Spline Curve');
legend;
grid on;
xlabel('X');
ylabel('Y');
title('B-Spline Curve');
hold off;

% De Boor algorithm for B-spline curve evaluation
function point = de_boor(ctrl_pts, knots, d, t)
    % De Boor's algorithm to evaluate B-spline curve at parameter t
    % ctrl_pts - control points [Nx2]
    % knots - knot vector
    % d - degree of the B-spline
    % t - parameter value to evaluate
    n = size(ctrl_pts, 1); % Number of control points
    i = find(knots <= t, 1, 'last'); % Knot span
    if isempty(i) || i < d+1
        i = d+1; % Ensure valid index for i
    end

    % Initialize de Boor points
    points = ctrl_pts((i-d):i, :);

    % Iterative de Boor computation
    for r = 1:d
        for j = (d+1):-1:(r+1)
            k_idx = i - d + j - 1; % Knot index for alpha computation
            alpha = (t - knots(k_idx)) / (knots(k_idx + d - r + 1) - knots(k_idx));
            points(j, :) = (1 - alpha) * points(j-1, :) + alpha * points(j, :);
        end
    end

    % Output the evaluated point
    point = points(end, :);
end
