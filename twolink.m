close all
% Parameters
l1 = 1; l2 = 1; % Link lengths
m1 = 1; m2 = 1; % Masses
g = 0;

% Trajectory definition (e.g., sinusoidal)
t = linspace(0, 5, 501); % Time vector (5 seconds, 501 points)
K = -1;

% Desired Joint Trajectory
q_Des = [ones(1, length(t)) * pi/4; ones(1, length(t)) * pi/4];
 
% Initial joint angles and velocities
q_Act = [0; 0]; % Start at the first point
qd_Act = [0; 0]; % Initial velocity

% Controller gains
Kp = [1, 0; 0, 1]; % Proportional gain (adjust as needed)
Kd = [.1, 0; 0, .1]; % Derivative gain (adjust as needed)

% Simulation parameters
dt = t(2) - t(1); % Time step

% Preallocate memory
n_points = length(t);
trajectory_qAct = zeros(2, n_points); % Store joint positions for plotting
tau = zeros(n_points, 2); % Store torque
Force = zeros(n_points, 2); % Store force
Pos = zeros(n_points, 2); % Store positions

% Loop through the trajectory
for i = 1:n_points
    % Current desired position and velocity
    qd_Des = [0; 0]; % Assuming a constant velocity profile
    
    % Error in position and velocity
    e = q_Des(:, i) - q_Act;
    eDot = qd_Des - qd_Act;
    
    % Compute dynamics
    M = mass_matrix(q_Act(1), q_Act(2), l1, l2, m1, m2);
    G = gravity_vector(q_Act(1), q_Act(2), l1, l2, m1, m2, 0, g);
    C = coriolis_matrix(q_Act(1), q_Act(2), qd_Act(1), qd_Act(2), l1, l2, m1, m2);
    
    % PD control with gravity compensation
    Torque = Kp * e * exp(-e'*Kp*e) + Kd * eDot;
    
    % Update joint dynamics using Euler integration
    qd_Act = qd_Act + dt * (M \ (Torque - C)); % Update velocity
    q_Act = q_Act + dt * qd_Act; % Update position
    
    % Store for visualization
    J_current = jacobian_2link(q_Act(1), q_Act(2), l1, l2);
    F = (-K * (inv(J_current') * Torque))';
    tau(i, :) = Torque'; % Save torque
    Force(i, :) = F; % Save force
    Pos(i, :) = forward_kinematics(q_Act(1), q_Act(2), l1, l2); % Save position
    trajectory_qAct(:, i) = q_Act; % Save joint angles
end
% Plot trajectory
figure(1); hold on; grid on;
plot(t, q_Des(1, :), '--', 'DisplayName', 'q1 Desired');
plot(t, q_Des(2, :), '--', 'DisplayName', 'q2 Desired');
plot(t, trajectory_qAct(1, :), '-', 'DisplayName', 'q1 Actual');
plot(t, trajectory_qAct(2, :), '-', 'DisplayName', 'q2 Actual');
xlabel('Time (s)');
ylabel('Joint Angles (rad)');
title('Joint Space Trajectory Following');
legend show;

figure(2); hold on; grid on;
plot(q_Des(1, 1), q_Des(2, 1), '*', 'DisplayName', 'Equilibrium Point');
plot(trajectory_qAct(1, 1), trajectory_qAct(2, 1), '*', 'DisplayName', 'Starting Point');
quiver(trajectory_qAct(1, :), trajectory_qAct(2, :), tau(:, 1)', tau(:, 2)');
legend show;

% Cartesian Space
figure(3); hold on; grid on;
plot(t, Pos(:,1), '--', 'DisplayName', 'q1 Desired');
plot(t, Pos(:,1), '--', 'DisplayName', 'q2 Desired');
xlabel('Time (s)');
ylabel('Joint Angles (rad)');
title('Joint Space Trajectory Following');
legend show;


figure(4); hold on; grid on;
plot(Pos(1,1),Pos(1,2),'*', 'DisplayName', 'First Position')
plot(Pos(end,1),Pos(end,2),'*', 'DisplayName', 'Last Position')
quiver(Pos(:,1), Pos(:,2), Force(:, 1), Force(:, 2));
legend show;

