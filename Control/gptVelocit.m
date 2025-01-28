clear; clc;

% Define desired trajectory with velocities
qDes = [ -0.4986,  2.5681, 0.1, 0.2;   % [q1, q2, q1_dot, q2_dot]
          0.5371,  1.5108, 0.05, -0.1 ];

% Optimization setup
initParams = [4 8  2 20 45]; % Initial guess for [time, wn, bj, kj]

[init_T, init_Y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, initParams(3), initParams(1:2), qDes, initParams(4), initParams(5)), [0 initParams(2)], zeros(8, 1));

% Lower and upper boundaries 
lb = [1 2  1  1      2  ];   % Lower bounds
ub = [5 5  40   100  200 ];  % Upper bounds

% Objective Function
objectiveFunc = @(params) objectiveFunction(params, qDes);

% Run optimization
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'MaxIter', 200);
optimalParams = fmincon(objectiveFunc, initParams, [], [], [], [], lb, ub, [], options);

% Simulate with optimal parameters and plot results
[t, y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, optimalParams(3), optimalParams(1:2), qDes, optimalParams(4), optimalParams(5)), [0 optimalParams(2)], zeros(8, 1));
xAct = forward_kinematics(y(:, 5), y(:, 6), 1, 1);
xDes = forward_kinematics(qDes(:, 1), qDes(:, 2), 1, 1);
xInit = forward_kinematics(init_Y(:, 5), init_Y(:, 6), 1, 1);

figure(1); hold on;
plot(xInit(:, 1), xInit(:, 2), '-');
plot(xAct(:, 1), xAct(:, 2), '-');
plot(xDes(:, 1), xDes(:, 2), 'o-');
legend('Initial','Optimised', 'Desired');
title('Optimized Trajectory Tracking');
disp(['Optimized Parameters :', num2str(optimalParams)])

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

    % weights
    w1 = 2000;
    w2 = 1;

    % Calculate the error metric (e.g., sum of squared errors)
    distto1 = w1 * sum((y(:, 5:6) - qDes(1, 1:2)).^2,2) + w2 * sum(abs(time(1) - t)); 
    distto2 = w1 * sum((y(:, 5:6) - qDes(2, 1:2)).^2,2) + w2 * sum(abs(time(2) - t)); 
    error   = min(distto1) + min(distto2);
end

% myTwolinkwithprefilter function
function dxdt = myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj)
    zeta = 1.0;
    A = [zeros([2 2]) eye(2); 
         -eye(2)*wn^2 -eye(2)*2*zeta*wn];
    B = [zeros(2, 2); wn^2*eye(2)];

    % Desired position and velocity
    if t < time(1)
        qDes_t = qDes(1, 1:2)';
        qDes_dot_t = qDes(1, 3:4)';
    else
        qDes_t = qDes(2, 1:2)';
        qDes_dot_t = qDes(2, 3:4)';
    end
    
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

    % Prefilter dynamics including desired velocities
    dotx = A*x(1:4) + B*[qDes_t; qDes_dot_t];
    dxdt = [dotx; qd; qdd];
end  
