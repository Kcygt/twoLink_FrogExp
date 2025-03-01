% Optimization of three-link robot arm tracking
clear; clc;

% Define desired trajectory and Middle Points
qDes = [ -0.5  2.5  1.0;
          0.5  1.5  0.8];

qMid = [inverse_kinematics(0.4, 0.6, 1, 1, 1), ...
        inverse_kinematics(0.4, 0.7, 1, 1, 1), ...
        inverse_kinematics(0.4, 0.8, 1, 1, 1), ...
        inverse_kinematics(0.4, 0.9, 1, 1, 1), ...
        inverse_kinematics(0.4, 1.0, 1, 1, 1), ...
        inverse_kinematics(0.4, 1.1, 1, 1, 1), ...
        inverse_kinematics(0.4, 1.2, 1, 1, 1)];

%  Parameters
time = [10 20];       % time
wn = [.5 2];          % Prefilter Omega     
kj = [40 25 30];      % Spring  [q1 q2 q3]
bj = [10 30 20];      % Damping [q1 q2 q3]
wt = [400, 1, 1800];  % weights [qDes, Time, qMid]

% Optimization setup
initParams = [time wn bj kj]; % Initial guess for [time, wn, bj, kj]
% Run initial simulation with given parameters
initState = zeros(12, 1); % 12 states for a three-link robot

[init_T, init_Y] = ode45(@(t, x) myThreelinkwithprefilter(t, x, wn, initParams(1:2), qDes, bj, kj), ...
                        [0 initParams(2)], initState);



% Objective function
function error = objectiveFunction(params, qDes, wt, qMid)
    % Initial conditions
    x0 = zeros(12, 1);
    x0(1:3) = [qDes(1, 1); qDes(1, 2); qDes(1, 3)];
    
    % Simulate the system
    [t, y] = ode45(@(t, x) myThreelinkwithprefilter(t, x, params(3:5), params(1:2), qDes, params(6:8), params(9:11)), [0 params(2)], x0);

    % Compute error metrics 
    error = sum(wt .* sum((y(:, 7:9) - qDes).^2, 2));
end

% Dynamics function
function dxdt = myThreelinkwithprefilter(t, x, wn, time, qDes, bj, kj)
    zeta = 1;
    A1 = blkdiag([zeros(3), eye(3); -eye(3)*wn(1)^2, -eye(3)*2*zeta*wn(1)]);
    B1 = [zeros(3); wn(1)^2*eye(3)];
    
    A2 = blkdiag([zeros(3), eye(3); -eye(3)*wn(2)^2, -eye(3)*2*zeta*wn(2)]);
    B2 = [zeros(3); wn(2)^2*eye(3)];

    % Positions and velocities
    q = x(7:9);
    qd = x(10:12);

    % Mass and inertia
    L = [1, 1, 1];
    M = diag([sum(L), sum(L), sum(L)]); % Approximate mass matrix
    V = -diag(bj) * qd - diag(kj) * (q - x(1:3)); 

    % Compute qdd
    qdd = M \ V;
    
    % Control logic
    if t < time(1)
        dotx = A1*x(1:6) + B1*qDes(1, :)';
    else
        dotx = A2*x(1:6) + B2*qDes(2, :)';
    end
    
    dxdt = [dotx; qd; qdd];
end

function qDes = inverse_kinematics(x, y, l1, l2, l3)
    r = sqrt(x.^2 + y.^2);
    
    % Ensure cos_q3 is within valid range
    cos_q3 = (r.^2 - (l1^2 + l2^2 + l3^2)) / (2 * l2 * l3);
    cos_q3 = min(max(cos_q3, -1), 1);  % Clamping to [-1,1]

    sin_q3 = sqrt(1 - cos_q3.^2); 
    q3 = atan2(sin_q3, cos_q3);
    
    % Compute q2
    k1 = l1 + l2*cos_q3 + l3*cos(q3);
    k2 = l2*sin_q3 + l3*sin(q3);
    q2 = atan2(k2, k1);

    % Compute q1
    phi = atan2(y, x);
    q1 = phi - q2;
    
    qDes = [q1; q2; q3];
end


% Forward kinematics
function P = forward_kinematics(q1, q2, q3, l1, l2, l3)
    x = l1*cos(q1) + l2*cos(q1 + q2) + l3*cos(q1 + q2 + q3);
    y = l1*sin(q1) + l2*sin(q1 + q2) + l3*sin(q1 + q2 + q3);
    P = [x, y];
end
