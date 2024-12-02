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

% Controller gains
K = 100; % Proportional gain
B = 40; % Derivative gain

% Initial joint angles and velocities
q0 = [0; 0]; % Initial joint angles
qd0 = [0; 0]; % Initial joint velocities
initial_state = [q0; qd0];


% Solve the ODE
[t, state] = ode45(@(t, x) robot_dynamics(t, x, l1, l2, m1, m2, g, K, B), t, initial_state);

% Extract results
qPos = state(:, 1:2); % Joint positions
qVel = state(:, 3:4); % Joint velocities

% Compute end-effector positions
xPos = zeros(length(t), 2);
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
    xPos(i, :) = forward_kinematics(q(1), q(2), l1, l2);
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
plot(t, xPos(:, 1), '--', 'DisplayName', 'End Effector Desired X');
plot(t, xPos(:, 2), '--', 'DisplayName', 'End Effector Desired Y');
xlabel('Time (s)');
ylabel('End Effector Position');
title('Cartesian Space Trajectory Following');
legend show;

figure(4); hold on; grid on;
step = 50:501;
plot(xPos(1, 1), xPos(1, 2), '*', 'DisplayName', 'First Position');
plot(xPos(end, 1), xPos(end, 2), '*', 'DisplayName', 'Last Position');
quiver(xPos(step, 1), xPos(step, 2), Force(step, 1), Force(step, 2));
legend show;

