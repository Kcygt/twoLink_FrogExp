clear; clc;

% Define desired trajectory and Middle Points
qDes = [0.1914, -0.0445, 0.3336];  % Desired final position

qMid = zeros(1,3); 
% qMid(1,:) = IK(0.02, 0, 0.01);
qMid(1,:) = IK(0.05, 0, 0.0);
% qMid(3,:) = IK(0.04, 0, 0.03);

% Parameters (fixed)
kj = [30 20 40];  % Spring constants
bj = [5 3 4];  % Damping constants
wt = [1000, 0.01, 2000000];  % Weights [qDes, Time, qMid]

% Optimization setup (only optimizing time and wn)
initParams = [10,5.5991   ,   2.5611,      7.5633];  % Initial guess for [time, wn(1), wn(2), wn(3)]

% Initial simulation with fixed parameters
[init_T, init_Y] = ode15s(@(t, x) myTwolinkwithprefilter(t, x, initParams(2:4), initParams(1), qDes, bj, kj), [0 initParams(1)], zeros(12, 1));
[xInit, yInit, zInit] = FK(init_Y(:,7), init_Y(:,8), init_Y(:,9));
[xD, yD, zD] = FK(init_Y(:,1), init_Y(:,2), init_Y(:,3));

% Define lower and upper boundaries for optimization
lb = [0, 1, 1, 1];  % Lower bounds for [time, wn(1), wn(2), wn(3)]
ub = [4, 10, 10, 10];  % Upper bounds for [time, wn(1), wn(2), wn(3)]

% Objective Function
objectiveFunc = @(params) objectiveFunction(params, qDes, wt, qMid, bj, kj);

% Run optimization
options = optimset('Display', 'iter', 'TolFun', 1e-7, 'MaxIter', 200);
optimalParams = fmincon(objectiveFunc, initParams, [], [], [], [], lb, ub, [], options);

% Simulate with optimal parameters
[t, y] = ode15s(@(t, x) myTwolinkwithprefilter(t, x, optimalParams(2:4), optimalParams(1), qDes, bj, kj), [0 optimalParams(1)], zeros(12, 1));

% Output
[xAct, yAct, zAct] = FK(y(:,7), y(:,8), y(:,9));
[xDes, yDes, zDes] = FK(qDes(1), qDes(2), qDes(3));

% Plot results
figure(1); hold on; grid on;
plot(xInit, zInit, '-.');
plot(xAct, zAct, '-');
plot(xDes, zDes, '*');
plot(0.02, 0.01, '*');
plot(0.03, 0.02, '*');
plot(0.04, 0.03, '*');

xlabel('X axis'); ylabel('Z axis');
legend('Initial', 'Optimized', 'Desired');
title('Initial vs Optimized Trajectory Tracking');

disp(['Optimized Parameters: ', num2str(optimalParams)]);

% Objective function
function error = objectiveFunction(params, qDes, wt, qMid, bj, kj)
    % Initial state vector
    x0 = zeros(12, 1);
    x0(1:3) = qDes;  % Set initial joint positions to qDes
    
    % Simulate the system with current parameters
    [t, y] = ode15s(@(t, x) myTwolinkwithprefilter(t, x, params(2:4), params(1), qDes, bj, kj), [0 params(1)], x0);

    % Calculate the distance to the desired trajectory
    distTo1 = min(sum((y(:, 7:9) - qDes).^2, 2) + sum((params(1) - t).^2, 2)); 

    % Calculate the distance to the middle points
    distMid = sum(arrayfun(@(i) min(sum((y(:, 7:9) - qMid(i, :)).^2, 2)), 1:size(qMid,1)));

    % Objective error
    error = wt(1) * distTo1 + wt(2) * params(1) + wt(3) * distMid;
end

% myTwolinkwithprefilter function
function dxdt = myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj)
    zeta = 1;  % Damping ratio
    A = [zeros(3,3), eye(3);
         -eye(3)*diag(wn).^2, -eye(3)*2*zeta*diag(wn)];
    B = [zeros(3,3); diag(wn).^2];

    q   = x(7:9);  % Joint positions
    qd  = x(10:12);  % Joint velocities

    Kp = -diag(kj);  % Spring constants (proportional gain)
    Kd = -diag(bj);  % Damping constants (derivative gain)

    % Compute mass, Coriolis, and gravity matrices
    [M, C, G] = compute_M_C_G(q(1), q(2), q(3), qd(1), qd(2), qd(3));

    % Joint accelerations
    qdd = M \ (C * qd + Kd * (qd - x(4:6)) + Kp * (q - x(1:3)));

    % Differential equations
    dxdt = [A*x(1:6) + B*qDes(:); qd; qdd];
end

% Forward Kinematics function
function [x, y, z] = FK(q1, q2, q3)
    l1 = 0.208;  % Length of first link
    l2 = 0.168;  % Length of second link
    x = sin(q1) .* (l1 * cos(q2) + l2 * sin(q3));
    y = l2 - l2 * cos(q3) + l1 * sin(q2);
    z = -l1 + cos(q1) .* (l1 * cos(q2) + l2 * sin(q3));
end

% Inverse Kinematics function
function [q1, q2, q3] = IK(x, y, z)
    l1 = 0.208;  % Length of first link
    l2 = 0.168;  % Length of second link
    q1 = atan2(x, z + l1);

    R = sqrt(x^2 + (z + l1)^2);
    r = sqrt(x^2 + (y - l2)^2 + (z + l1)^2);
    Beta = atan2(y - l2, R);
    Gamma = acos((l1^2 + r^2 - l2^2) / (2 * l1 * r));
    q2 = Gamma + Beta;

    Alpha = acos((l1^2 + l2^2 - r^2) / (2 * l1 * l2));
    q3 = q2 + Alpha - pi/2;
end
