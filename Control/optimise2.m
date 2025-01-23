% Optimization of two-link robot arm tracking
clear; clc;

% Define desired trajectory
[fX, fY, cartesianX, cartesianY] = defineSquare(0.9, 0.9, 1);
qDes = inverse_kinematics(cartesianX, cartesianY, 1, 1)';  % Desired joint angles

% Optimization setup
initParams = [3, 6, 9, 12, 15,     1   20,   45]; % Initial guess for [time, wn, bj, kj]
% initialParams =   [  1.9142    3.3451    4.0022    6.0404    8.0008    4.8515   10.4579   47.2359];
initParams = [    1.7300    3.9999    5.0020    6.0017    8.0005    4.5564   12.3538   49.7463
]; % Initial guess for [time, wn, bj, kj]

[init_T, init_Y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, initParams(6), initParams(1:5), qDes, initParams(7), initParams(8)), [0 initParams(5)], zeros(8, 1));

% Lower and upper boundaries 
tStep = 2;
lb = [0,     tStep,   2*tStep, 3*tStep, 4*tStep,   0.1,    1,  10  ];   % Lower bounds
ub = [tStep, 2*tStep, 3*tStep, 4*tStep, 5*tStep,    5,   50,  100 ];     % Upper bounds

% Objective Function
objectiveFunc = @(params) objectiveFunction(params, qDes);

% Run optimization
options = optimset('Display', 'iter', 'TolFun', 1e-7, 'MaxIter', 200);
optimalParams = fmincon(objectiveFunc, initParams, [], [], [], [], lb, ub, [], options);

% Simulate with optimal parameters and plot results
[t, y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, optimalParams(6), optimalParams(1:5), qDes, optimalParams(7), optimalParams(8)), [0 optimalParams(5)], zeros(8, 1));
xAct = forward_kinematics(y(:, 5), y(:, 6), 1, 1);
xDes = forward_kinematics(qDes(:, 1), qDes(:, 2), 1, 1);
xInit = forward_kinematics(init_Y(:, 5), init_Y(:, 6), 1, 1);

figure(1); hold on;
plot(xInit(:, 1), xInit(:, 2), '-');
plot(xAct(:, 1), xAct(:, 2), '-');
plot(xDes(:, 1), xDes(:, 2), 'o-');
legend('Initial','Optimised', 'Desired');
title('Optimized Trajectory Tracking');


% Objective function
function error = objectiveFunction(params, qDes)
    time = [params(1), params(2), params(3), params(4), params(5)];
    
    wn = params(6);
    bj = params(7);
    kj = params(8);
   
    % Initial conditions
    x0 = zeros(8, 1);
    x0(1:2) = [qDes(1, 1); qDes(1, 2)];
    
    % Simulate the system
    [t, y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj), [0 time(end)], x0);

    % weights
    w1 = 400;
    w2 = 1;

    % Calculate the error metric (e.g., sum of squared errors)
    distto1 = w1 * sum((y(:,5:6) - qDes(1,:)).^2,2) + w2 * sum(abs(time(1) - t)); 
    distto2 = w1 * sum((y(:,5:6) - qDes(2,:)).^2,2) + w2 * sum(abs(time(2) - t));
    distto3 = w1 * sum((y(:,5:6) - qDes(3,:)).^2,2) + w2 * sum(abs(time(3) - t));
    distto4 = w1 * sum((y(:,5:6) - qDes(4,:)).^2,2) + w2 * sum(abs(time(4) - t));
    distto5 = w1 * sum((y(:,5:6) - qDes(5,:)).^2,2) + w2 * sum(abs(time(5) - t));
    error   = min(distto1) + min(distto2) + min(distto3) + min(distto4) + min(distto5);
    % error   =min(distto3);
end

% myTwolinkwithprefilter function
function dxdt = myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj)
    zeta = 1.0;
    A = [zeros([2 2]) eye(2); -eye(2)*wn^2 -eye(2)*2*zeta*wn];
    B = [0 0; 0 0; wn^2 0; 0 wn^2];
    
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
    elseif t < time(2)
        dotx = A*x(1:4) + B*qDes(2, :)';
    elseif t < time(3)
        dotx = A*x(1:4) + B*qDes(3, :)';
    elseif t < time(4)
        dotx = A*x(1:4) + B*qDes(4, :)';
    else
        dotx = A*x(1:4) + B*qDes(5, :)';
    end
    
    dxdt = [dotx; qd; qdd];
end

function [sX, sY,cX,cY] = defineSquare(x_c, y_c, a)

    % Calculate half side length
    half_side = a / 2;

    % Define vertices
    cX = [ x_c - half_side,  x_c - half_side,  x_c + half_side, x_c + half_side, x_c - half_side];
    cY = [ y_c - half_side,  y_c + half_side,  y_c + half_side, y_c - half_side, y_c - half_side];

    mX = [(cX(2) + cX(1))/2, (cX(3) + cX(2))/2, (cX(4) + cX(3))/2, (cX(5) + cX(4))/2];
    mY = [(cY(2) + cY(1))/2, (cY(3) + cY(2))/2, (cY(4) + cY(3))/2, (cY(5) + cY(4))/2];
    
    sX = [cX(1), mX(1),cX(2), mX(2),cX(3), mX(3),cX(4), mX(4),cX(5)];
    sY = [cY(1), mY(1),cY(2), mY(2),cY(3), mY(3),cY(4), mY(4),cY(5)];
    
end
