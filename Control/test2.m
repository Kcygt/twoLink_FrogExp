% Optimization of two-link robot arm tracking
clear; clc;

% Define desired trajectory
qDes = [ -0.4986    2.5681;
          0.5371    1.5108 ];

% Optimization setup
initParams = [10 20   1 20 45]; % Initial guess for [time, wn, bj, kj]

[init_T, init_Y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, initParams(3), initParams(1:2), qDes, initParams(4), initParams(5)), [0 initParams(2)], zeros(8, 1));

% Lower and upper boundaries 
lb = [0 0    1.5   1    2  ];   % Lower bounds
ub = [3 6   50  200  500 ];  % Upper bounds

% Objective Function
objectiveFunc = @(params) objectiveFunction(params, qDes);

% Run optimization
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'MaxIter', 200);
optimalParams = fmincon(objectiveFunc, initParams, [], [], [], [], lb, ub, [], options);

% Simulate with optimal parameters and plot results
[t, y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, optimalParams(3), optimalParams(1:2), qDes, optimalParams(4), optimalParams(5)), ...
               [0 optimalParams(2)], zeros(8, 1));

% Compute the actual trajectory
xAct = forward_kinematics(y(:, 5), y(:, 6), 1, 1);

% Compute the desired straight-line trajectory for visualization
num_points = size(y, 1);
xDesLine = [linspace(qDes(1,1), qDes(2,1), num_points)', linspace(qDes(1,2), qDes(2,2), num_points)'];

% Plot results
figure(1); hold on;
plot(xAct(:, 1), xAct(:, 2), 'b-', 'LineWidth', 2); % Optimized trajectory
plot(xDesLine(:, 1), xDesLine(:, 2), 'r--', 'LineWidth', 2); % Ideal straight-line trajectory
plot(qDes(:, 1), qDes(:, 2), 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Desired waypoints
legend('Optimized Path', 'Ideal Straight Line', 'Waypoints');
title('Optimized Trajectory Tracking');
grid on;

% Objective function
function error = objectiveFunction(params, qDes)
    time = [params(1) params(2)];
    wn = params(3);
    bj = params(4);
    kj = params(5);
    
    % Initial conditions
    x0 = zeros(8, 1);
    x0(1:2) = [qDes(1, 1); qDes(1, 2)];
    
    % Simulate the system
    [t, y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj), [0 time(end)], x0);
    
    % Extract the actual trajectory
    q_actual = y(:, 5:6);
    
    % Define the straight line equation
    q1_start = qDes(1, 1);
    q2_start = qDes(1, 2);
    q1_end = qDes(2, 1);
    q2_end = qDes(2, 2);
    
    % Line equation coefficients (Ax + By + C = 0)
    A = q2_start - q2_end;
    B = q1_end - q1_start;
    C = q1_start * q2_end - q1_end * q2_start;
    
    % Distance of each actual point from the line
    line_distances = abs(A * q_actual(:, 1) + B * q_actual(:, 2) + C) ./ sqrt(A^2 + B^2);
    
    % Weights for different penalties
    w1 = 1000;  % Weight for deviation at start and end points
    w2 = 500;   % Weight for deviation from the straight line
    
    % Error metric: minimize deviation from the line and ensure reaching desired points
    error = w1 * sum((q_actual(1, :) - qDes(1, :)).^2) + ...
            w1 * sum((q_actual(end, :) - qDes(2, :)).^2) + ...
            w2 * sum(line_distances.^2);
end


% myTwolinkwithprefilter function
function dxdt = myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj)
    zeta = 1.0;
    A = [zeros([2 2]) eye(2); -eye(2)*wn^2 -eye(2)*2*zeta*wn];
    B = [0 0; 0 0; wn^2 0; 0 wn^2];
    
    % Actual position and velocity
    q = x(5:6);
    qd = x(7:8);
    q1p = x(7); q2p = x(8);
    q1 = x(5); q2 = x(6);
    
    % Robot constants
    L_1 = 1; L_2 = 1; m_1 = 1; m_2 = 1;
    ka = L_2^2 * m_2;
    kb = L_1 * L_2 * m_2;
    kc = L_1^2 * (m_1 + m_2);
    
    M = [ka + 2*kb*cos(q2) + kc, ka + kb*cos(q2);
         ka + kb*cos(q2), ka];
    V = ka*sin(q2)*([0 -1; 1 0] * [q1p^2; q2p^2] + [-2*q1p*q2p; 0]);
    
    Numerator = V + [-bj 0; 0 -bj]*qd + [-kj 0; 0 -kj]*(q - x(1:2));
    qdd = M\Numerator;
    if t < time(1)
        dotx = A*x(1:4) + B*qDes(1, :)';
    else
      dotx = A*x(1:4) + B*qDes(2, :)';
    end
    dxdt = [dotx; qd; qdd];
end
