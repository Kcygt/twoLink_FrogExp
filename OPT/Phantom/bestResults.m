% Optimization of two-link robot arm tracking
clear; clc;

% Define desired trajectory and Middle Points
qDes = [ -0.4986    2.5681;
          0.5371    1.5108 ];

qMid = [inverse_kinematics(0.4, 0.6, 1, 1), ...
        inverse_kinematics(0.4, 0.7, 1, 1), ...
        inverse_kinematics(0.4, 0.8, 1, 1), ...
        inverse_kinematics(0.4, 0.9, 1, 1), ...
        inverse_kinematics(0.4, 1.0, 1, 1), ...
        inverse_kinematics(0.4, 1.1, 1, 1), ...
        inverse_kinematics(0.4, 1.2, 1, 1)];

%  Parameters
time = [10 20];      % time
wn = [20 15];        % Prefilter Omega     
kj = [40 25];        % Spring  [q1 q2]
bj = [10 30];        % Damping [q1 q2]
wt = [400, 1, 1800];  % weights [qDes, Time, qMid]

% Optimization setup
initParams = [time  wn bj kj]; % Initial guess for [time, wn, bj, kj]

[init_T, init_Y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, wn, initParams(1:2), qDes, bj, kj),   [0 initParams(2)], zeros(8, 1));

% Lower and upper boundaries 
lb = [0 0   1.5 1.5   10  10    2   2     ];   % Lower bounds
ub = [2 6   10   10     200 200   200 200 ];     % Upper bounds

% Objective Function
objectiveFunc = @(params) objectiveFunction(params, qDes, wt, qMid);

% Run optimization
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'MaxIter', 400);
optimalParams = fmincon(objectiveFunc, initParams, [], [], [], [], lb, ub, [], options);

% Simulate with optimal parameters and plot results
[t, y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, optimalParams(3:4), optimalParams(1:2), qDes, optimalParams(5:6), optimalParams(7:8)), [0 optimalParams(2)], zeros(8, 1));

% Output
xAct = forward_kinematics(y(:, 5), y(:, 6), 1, 1);
xDes = forward_kinematics(qDes(:, 1), qDes(:, 2), 1, 1);
xInit = forward_kinematics(init_Y(:, 5), init_Y(:, 6), 1, 1);

% Plotting
% Desired, Actual and Optimised Data
figure(1); hold on; grid on;
plot(xInit(:, 1), xInit(:, 2), '-');
plot(xAct(:, 1), xAct(:, 2), '-');
plot(xDes(:, 1), xDes(:, 2), 'o-');
plot(0.4,0.6, '*',0.4,0.8, '*',0.4,0.9, '*',0.4,1.2, '*'); 
xlabel('X axis'); ylabel('Y axis');
legend('Initial','Optimised', 'Desired');
title('Optimized Trajectory Tracking');
disp(['Optimized Parameters :', num2str(optimalParams)])

% Objective function
function error = objectiveFunction(params, qDes,wt,qMid)
    
    % Initial conditions
    x0 = zeros(8, 1);
    x0(1:2) = [qDes(1, 1); qDes(1, 2)];
    
    % Simulate the system
    [t, y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, params(3:4), params(1:2), qDes, params(5:6), params(7:8)),   [0 params(2)], x0);

    % Calculate the error metric 
    distto1 = min(sum((y(:, 5:6) - qDes(1,:)).^2,2) + sum((params(1) - t).^2,2) ); 
    distto2 = min(sum((y(:, 5:6) - qDes(2,:)).^2,2) + sum((params(2) - t).^2,2) ); 

    distMid1 = min(sum((y(:, 5:6) - qMid(:,1)').^2,2));  
    distMid2 = min(sum((y(:, 5:6) - qMid(:,2)').^2,2));       
    distMid3 = min(sum((y(:, 5:6) - qMid(:,3)').^2,2));       
    distMid4 = min(sum((y(:, 5:6) - qMid(:,4)').^2,2));      
    distMid5 = min(sum((y(:, 5:6) - qMid(:,5)').^2,2));      
    distMid6 = min(sum((y(:, 5:6) - qMid(:,6)').^2,2));      
    distMid7 = min(sum((y(:, 5:6) - qMid(:,7)').^2,2));      
    
    % time1 = min(sum((params(1) - t).^2,2));
    % time2 = min(sum((params(2) - t).^2,2));
    time1 = params(1);
    time2 = params(2);

    error   = wt(1) * distto1  + wt(1) * distto2  + ...  % Desired
              wt(2) * time1    + wt(2) * time2    + ...  % time
              wt(3) * distMid1 + wt(3) * distMid2 + ...  % Mid-point
              wt(3) * distMid3 + wt(3) * distMid4 + ...  % Mid-point
              wt(3) * distMid5 + wt(3) * distMid6 + ...  % Mid-point
              wt(3) * distMid7;                          % Mid-point


    % distto5 = 5000 * sum((y(:, 5:6) - qMid3'),2) + w2 * (sum(  (   (time(1) + (time(2) - time(1))/2 ) - t).^2   ,2));

end

% myTwolinkwithprefilter function
function dxdt = myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj)
    zeta = 1;
    A1 = [zeros([2 2]) eye(2); -eye(2)*wn(1)^2 -eye(2)*2*zeta*wn(1)];
    B1 = [0 0; 0 0; wn(1)^2 0; 0 wn(1)^2];
    
    A2 = [zeros([2 2]) eye(2); -eye(2)*wn(2)^2 -eye(2)*2*zeta*wn(2)];
    B2 = [0 0; 0 0; wn(2)^2 0; 0 wn(2)^2];
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
    
    Numerator = V + [-bj(1) 0; 0 -bj(2)]*qd + [-kj(1) 0; 0 -kj(2)]*(q - x(1:2));
    qdd = M\Numerator;
    if t < time(1)
        dotx = A1*x(1:4) + B1*qDes(1, :)';
    else
        dotx = A2*x(1:4) + B2*qDes(2, :)';
    end
    dxdt = [dotx; qd; qdd];
end

function qDes = inverse_kinematics(x, y, l1, l2)
    r = sqrt(x.^2 + y.^2);
    % if r > (l1 + l2) || r < abs(l1 - l2)
    %     error('Target point is out of reach');
    % end
    cos_q2 = (r.^2 - l1.^2 - l2.^2) / (2 * l1 * l2);
    sin_q2 = sqrt(1 - cos_q2.^2); 
    q2 = atan2(sin_q2, cos_q2);
    phi = atan2(y, x);
    psi = atan2(l2 * sin(q2), l1 + l2 * cos(q2));
    q1 = phi - psi;
    qDes = [q1; q2];
end

% Function Definitions
function P = forward_kinematics(q1, q2, l1, l2)
    x = l1*cos(q1) + l2*cos(q1 + q2);
    y = l1*sin(q1) + l2*sin(q1 + q2);
    P = [x, y];
end