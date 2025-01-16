% Define the optimization problem
objective_function = @(params) optimizeTrajectory(params);

% Set initial guesses for the parameters
initial_guess = [2, 40, 20, 1]; % [wn, kj, bj, total_time]

% Set bounds for the parameters (optional)
lb = [0.1, 10, 1, 0.1];   % lower bounds for [wn, kj, bj, total_time]
ub = [10, 100, 100, 5];   % upper bounds for [wn, kj, bj, total_time]

% Use fminunc for unconstrained optimization
options = optimset('Display', 'iter');
optimal_params = fminunc(objective_function, initial_guess, options);

% Display the optimized parameters
disp('Optimized parameters:');
disp(optimal_params);

% Define the function for optimizing the trajectory
function error = optimizeTrajectory(params)
    % Extract parameters
    wn = params(1);        % Natural frequency
    kj = params(2);        % Controller gain
    bj = params(3);        % Damping coefficient
    total_time = params(4); % Total time for simulation
    
    % Define square parameters (adjust if necessary)
    [fX, fY, cartesianX, cartesianY] = defineSquare(0.9, 0.9, 1);
    
    % Inverse kinematics to find initial joint angles
    qDes = inverse_kinematics(cartesianX, cartesianY, 1, 1)';
    
    % Set initial conditions for simulation
    x0 = zeros(8, 1);
    x0(1:2) = [qDes(1,1); qDes(1,2)];
    
    % Create time intervals for different phases
    time = linspace(0, total_time, 5);  % This creates 5 time points evenly spaced
    
    % Run the simulation with the given parameters
    [t, y] = ode113(@(t,x) myTwolinkwithprefilter(t, x, wn, time, qDes), [0 total_time], x0);

    % Forward kinematics to get the actual trajectory
    xAct = forward_kinematics(y(:, 5), y(:, 6), 1, 1);
    xDes = forward_kinematics(qDes(:, 1), qDes(:, 2), 1, 1);
    
    % Calculate error between actual and desired trajectories
    error = sum((xAct(:,1) - xDes(:,1)).^2 + (xAct(:,2) - xDes(:,2)).^2);
end

% Define the function to compute the dynamics of the system
function dxdt = myTwolinkwithprefilter(t, x, wn, time, qDes)
    zeta = 1.0;
    A = [zeros([2 2]) eye(2); -eye(2)*wn^2 -eye(2)*2*zeta*wn]; % note the wn^2 !!
    B = [0 0; 0 0; wn^2 0; 0 wn^2];
    
    % Set up a two-link arm
    q = x(5:6);
    qd = x(7:8);
    q1p = x(7); q2p = x(8);
    q1 = x(5); q2 = x(6);
    
    bj = 20; % energy dissipation term
    kj = 40; % controller gain
    L_1 = 1;
    L_2 = 1;
    m_1 = 1;
    m_2 = 1;
    
    % Derived constants
    ka = L_2^2 * m_2;
    kb = (L_1 * L_2 * m_2);
    kc = L_1^2 * (m_1 + m_2);
    
    M = [ka + 2*kb*cos(q2) + kc  ka + kb*cos(q2);
         ka + kb*cos(q2)        ka];
    V = ka * sin(q2) * ([0 -1; 1 0] * [q1p^2; q2p^2] + [-2*q1p*q2p; 0]);
    
    Numerator = V + [-bj 0; 0 -bj] * qd + [-kj 0; 0 -kj] * (q - x(1:2));
    qdd = M \ Numerator;
    
    % Define control input based on time
    if t < time(1)
        dotx = A * x(1:4) + B * qDes(1,:)';
    elseif t < time(2)
        dotx = A * x(1:4) + B * qDes(2,:)';
    elseif t < time(3)
        dotx = A * x(1:4) + B * qDes(3,:)';
    elseif t < time(4)
        dotx = A * x(1:4) + B * qDes(4,:)';
    else
        dotx = A * x(1:4) + B * qDes(5,:)';
    end
    
    dxdt = [dotx; qd; qdd];
end

% Define inverse kinematics function
function qDes = inverse_kinematics(x, y, l1, l2)
    r = sqrt(x.^2 + y.^2);
    cos_q2 = (r.^2 - l1.^2 - l2.^2) / (2 * l1 * l2);
    sin_q2 = sqrt(1 - cos_q2.^2); 
    q2 = atan2(sin_q2, cos_q2);
    phi = atan2(y, x);
    psi = atan2(l2 * sin(q2), l1 + l2 * cos(q2));
    q1 = phi - psi;
    qDes = [q1; q2];
end

% Define forward kinematics function
function x = forward_kinematics(q1, q2, l1, l2)
    x = [l1 * cos(q1) + l2 * cos(q1 + q2), l1 * sin(q1) + l2 * sin(q1 + q2)];
end

% Define function to define the square trajectory
function [sX, sY, cX, cY] = defineSquare(x_c, y_c, a)
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
