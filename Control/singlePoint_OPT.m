% Optimization of two-link robot arm tracking
clear; clc;

% Define desired trajectory
qDes = [ 0.5371    1.5108 ];

% Optimization setup
initParams = [20  1 20 45]; % Initial guess for [time, wn, bj, kj]

[init_T, init_Y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, initParams(2), initParams(1), qDes, initParams(3), initParams(4)), [0 initParams(1)], zeros(8, 1));

% Lower and upper boundaries 
lb = [0   1.5   1    2  ];   % Lower bounds
ub = [15  50  200  500 ];  % Upper bounds

% Objective Function

% Run optimization
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'MaxIter', 200);
optimalParams = fmincon(@(params) objectiveFunction(params, qDes), initParams, [], [], [], [], lb, ub, [], options);

% Simulate with optimal parameters and plot results
[t, y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, optimalParams(2), optimalParams(1), qDes, optimalParams(3), optimalParams(4)), [0 optimalParams(1)], zeros(8, 1));
xAct = forward_kinematics(y(:, 5), y(:, 6), 1, 1);
xDes = forward_kinematics(qDes(:, 1), qDes(:, 2), 1, 1);
xInit = forward_kinematics(init_Y(:, 5), init_Y(:, 6), 1, 1);

figure(1); hold on;
plot(xInit(:, 1), xInit(:, 2), '-');
plot(xAct(:, 1), xAct(:, 2), '-');
plot(xDes(:, 1), xDes(:, 2), 'o-');
plot(0.4,0.6, '*',0.4,0.8, '*',0.4,0.9, '*',0.4,1.2, '*'); 

legend('Initial','Optimised', 'Desired');
title('Optimized Trajectory Tracking');
disp(['Optimized Parameters :', num2str(optimalParams)])

% Objective function
function error = objectiveFunction(params, qDes)
    time =params(1);
    
    wn = params(2);
    bj = params(3);
    kj = params(4);
   
    % Initial conditions
    x0 = zeros(8, 1);
    x0(1:2) = [qDes(1, 1); qDes(1, 2)];
    
    % Simulate the system
    [t, y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj), [0 time], x0);

    % weights, could be done as a vector of weights
    w1 = 0;
    w2 = 0;
    w3 = 0;
    % w= [0.5, 1 , 5]; [qd_wt, time_wt, midpt_wt]

    % Mid Points
    qMid1 = inverse_kinematics(0.4, 0.6, 1, 1);
    qMid2 = inverse_kinematics(0.4, 0.8, 1, 1);
    qMid3 = inverse_kinematics(0.4, 0.9, 1, 1);
    qMid4 = inverse_kinematics(0.4, 1.2, 1, 1);

    % Calculate the error metric 
    distto1 = w1 * sum((y(:, 5:6) - qDes(1,:)).^2,2) + w2 * (sum((time(1) - t).^2,2)); 
    distto2 = 3000 * sum((y(:, 5:6) - qMid3').^2,2);
    error   = min(distto1) + min(distto2);

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
    dotx = A*x(1:4) + B*qDes(1, :)';
    dxdt = [dotx; qd; qdd];
end
