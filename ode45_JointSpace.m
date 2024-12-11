close all;
clear;

% Parameters
l1 = 1; l2 = 1; % Link lengths
m1 = 1; m2 = 1; % Masses
g = 0;

% Trajectory definition (e.g., sinusoidal)
t = linspace(0, 5, 501)'; % Time vector (5 seconds, 501 points)

% Desired Joint Trajectory
qDes = [ones(length(t), 1) * pi / 4, ones(length(t), 1) * pi / 4];
xDes = forward_kinematics(qDes(:,1), qDes(:,2), l1, l2);

% Controller gains
K = 100; % Proportional gain
B = 40; % Derivative gain

% Initial joint angles and velocities
q0 = [pi/3; -pi/3]; % Initial joint angles
qd0 = [0; 0]; % Initial joint velocities
initial_state = [q0; qd0];


% Solve the ODE
[t, state] = ode45(@(t, x) robot_dynamics(t, x, l1, l2, m1, m2, g, K, B), t, initial_state);

% Extract results
qPos = state(:, 1:2); % Joint positions
qVel = state(:, 3:4); % Joint velocities

% Compute end-effector positions
xAct = zeros(length(t), 2);
Force = zeros(length(t), 2);
tau = zeros(length(t), 2);

for i = 1:length(t)
    q = qPos(i, :)';
    qd = qVel(i, :)';
    J_current = jacobian_2link(q(1), q(2), l1, l2);

    % Recompute torque for visualization
    e = q - qDes(i, :)';
    eDot = qd - [0; 0];
    Torque = -K * e - B * eDot;

    % Store data
    xAct(i, :) = forward_kinematics(q(1), q(2), l1, l2);
    Force(i, :) = (-K * (inv(J_current') * Torque))';
    tau(i, :) = Torque';
end

% Plot trajectory
figure(1); hold on; grid on;
plot(t, qDes(:, 1), '--', 'DisplayName', 'q1 Desired');
plot(t, qDes(:, 2), '--', 'DisplayName', 'q2 Desired');
plot(t, qPos(:, 1), '-', 'DisplayName', 'q1 Actual');
plot(t, qPos(:, 2), '-', 'DisplayName', 'q2 Actual');
xlabel('Time (s)');
ylabel('Joint Angles (rad)');
title('Joint Space Trajectory Following');
legend show;

figure(2); hold on; grid on;
plot(qDes(1, 1), qDes(2, 1), '*', 'DisplayName', 'Equilibrium Point');
plot(qPos(1, 1), qPos(2, 1), '*', 'DisplayName', 'Starting Point');
quiver(qPos(:, 1), qPos(:, 2), tau(:, 1), tau(:, 2));
legend show;

% Cartesian Space
figure(3); hold on; grid on;
plot(xDes(:, 1),xDes(:, 2), 'o', 'DisplayName', 'End Effector Desired');
plot(xAct(:, 1),xAct(:, 2), '*', 'DisplayName', 'End Effector Actual');

xlabel('Time (s)');
ylabel('End Effector Position');
title('Cartesian Space Trajectory Following');
legend show;

figure(4); hold on; grid on;
step = 50:501;
plot(xAct(1, 1), xAct(1, 2), '*', 'DisplayName', 'First Position');
plot(xAct(end, 1), xAct(end, 2), '*', 'DisplayName', 'Last Position');
quiver(xAct(step, 1), xAct(step, 2), Force(step, 1), Force(step, 2));
legend show;


function dxdt = robot_dynamics(t, x, l1, l2, m1, m2, g, K, B)
    % Unpack state variables
    q = x(1:2);
    qd = x(3:4);

    % Desired trajectory
    q_des = [pi/4; pi/4]; % Fix for simplicity
    qd_des = [0; 0];

    % Error in position and velocity
    e = q - q_des;
    eDot = qd - qd_des;

    % Dynamics
    M = mass_matrix(q(1), q(2), l1, l2, m1, m2);
    G = gravity_vector(q(1), q(2), l1, l2, m1, m2, 0,0);
    C = coriolis_matrix(q(1), q(2), qd(1), qd(2), l1, l2, m1, m2);
    
    % PD control with gravity compensation
    Torque = -K * e - B * eDot;

    % State derivatives
    qdd = M \ (Torque - C * qd - G); % Joint accelerations
    dxdt = [qd; qdd];
end

% Function Definitions
function P = forward_kinematics(q1, q2, l1, l2)
    x = l1*cos(q1) + l2*cos(q1 + q2);
    y = l1*sin(q1) + l2*sin(q1 + q2);
    P = [x, y];
end

function qDes = inverse_kinematics(x, y, l1, l2)
    r = sqrt(x^2 + y^2);
    if r > (l1 + l2) || r < abs(l1 - l2)
        error('Target point is out of reach');
    end
    cos_q2 = (r^2 - l1^2 - l2^2) / (2 * l1 * l2);
    sin_q2 = sqrt(1 - cos_q2^2); 
    q2 = atan2(sin_q2, cos_q2);
    phi = atan2(y, x);
    psi = atan2(l2 * sin(q2), l1 + l2 * cos(q2));
    q1 = phi - psi;
    qDes = [q1, q2];
end

