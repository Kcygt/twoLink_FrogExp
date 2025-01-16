% Optimization of two-link robot arm tracking
clear; clc;

% Define constants and desired trajectory
[fX, fY, cartesianX, cartesianY] = defineSquare(0.9, 0.9, 1);
qDes = inverse_kinematics(cartesianX, cartesianY, 1, 1)';  % Desired joint angles
gain = 3.0;
% Optimization setup
timee = [ 2 4 6 8 10];
timeLow = [1, 2.1, 3.1, 4.1, 5.1];
timeUp  = [2,   3,   4,   5,   6];
initialParams = [timee,      3,   10,   20  ]; % Initial guess for [time, wn, bj, kj]
initialParams=[  1.1155    2.6454    3.9002    4.8978    5.8978    7.0857    9.9857   41.2005];
initialParams=[  1.1192    2.6448    3.9014    4.8941    5.8941    7.3920   17.4518   42.4887];
initialParams=[  1.3973    2.5916    3.6164    4.6204    5.6204    7.3605   17.4846   43.5925];


lb = [timeLow,                0.1,    1,    1  ];               % Lower bounds
ub = [timeUp,                 10,   40,  50 ];       % Upper bounds

% Use an anonymous function to pass qDes to the objective function
objectiveFunc = @(params) objectiveFunction(params, qDes);

% Run optimization
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'MaxIter', 2000);
optimalParams = fmincon(objectiveFunc, initialParams, [], [], [], [], lb, ub, [], options);

% Simulate with optimal parameters and plot results
[t, y] = ode113(@(t, x) myTwolinkwithprefilter(t, x, optimalParams(6), optimalParams(1:5), qDes, optimalParams(7), optimalParams(8)), [0 40], zeros(8, 1));
xAct = forward_kinematics(y(:, 5), y(:, 6), 1, 1);
xDes = forward_kinematics(qDes(:, 1), qDes(:, 2), 1, 1);

figure;
plot(xAct(:, 1), xAct(:, 2), '-');
hold on;
plot(xDes(:, 1), xDes(:, 2), 'o-');
legend('Actual', 'Desired');
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
    [t, y] = ode113(@(t, x) myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj), [0 30], x0);

    % index = zeros(length(time),1);
    % j = 1;
    % for i = 1:length(t)
    % 
    %     if abs(t(i) - time(j)) < 0.5
    %         index(j) = i;
    %         j = j + 1;
    %         if j == 6 
    %             break
    %         end
    %     end
    %     i = i+1;
    % 
    % end
    % index(end) = length(t)-1;
    
    % Compute the tracking error
    % xAct = forward_kinematics(y(:, 5), y(:, 6), 1, 1);
    % xDes = forward_kinematics(qDes(:, 1), qDes(:, 2), 1, 1);
    
    % Calculate the error metric (e.g., sum of squared errors)
    % error = abs(sum((y(:, 1) - qDes(1, 1)).^2 + (xAct(:, 2) - xDes(:, 2)).^2));

    distto1=sum((y(:,5:6)-qDes(1,:)).^2,2);
    distto2=sum((y(:,5:6)-qDes(2,:)).^2,2);
    distto3=sum((y(:,5:6)-qDes(3,:)).^2,2);
    distto4=sum((y(400:end,5:6)-qDes(4,:)).^2,2);
    error=min(distto1)+min(distto2)+min(distto3)+min(distto4);
    error=min(distto4);

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
