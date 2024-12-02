close all
clear
% Parameters
l1 = 1; l2 = 1; % Link lengths
m1 = 1; m2 = 1; % Masses
g = 0;

% Trajectory definition (e.g., sinusoidal)
t = linspace(0, 5, 501)'; % Time vector (5 seconds, 501 points)

% Desired Joint Trajectory
qDes = [ones(length(t),1) * pi/4, ones(length(t),1) * pi/4];
 
% Initial joint angles and velocities
qAct = [0 0]; % Start at the first point
qdAct = [0 0]; % Initial velocity


% Simulation parameters
dt = t(2) - t(1); % Time step

% Preallocate memory
n_points = length(t);
qPos = zeros(n_points,2); % Store joint positions for plotting
xPos = zeros(n_points, 2); % Store positions
qVel = zeros(n_points, 2);
xVel = zeros(n_points, 2);
tau = zeros(n_points, 2); % Store torque
Force = zeros(n_points, 2); % Store force


% Initialize qd_Des_Prev
qd_Des_Prev = [0, 0]; % Assuming initial desired velocity is zero
xd_Des_Prev = [0, 0];

% Co5troller gains
K = 50; % Proportional gain (adjust as needed)
B = 30; % Derivative gain (adjust as needed)

for i = 1:n_points
    % Update desired velocity
    qdDes = qd_Des_Prev - qdAct; 
    
    % Error in position and velocity
    e =  qAct - qDes(i,:);
    eDot = qdAct - qdDes ;
    
    % Compute dynamics
    M = mass_matrix(qAct(1), qAct(2), l1, l2, m1, m2);
    G = gravity_vector(qAct(1), qAct(2), l1, l2, m1, m2, 0, g);
    C = coriolis_matrix(qAct(1), qAct(2), qdAct(1), qdAct(2), l1, l2, m1, m2);
    
    % PD control with gravity compensation
    % Torque = K * -e * exp(-e'*0*e) + B * -eDot;
    Torque = K * -e  + B * -eDot;

    % Update joint dynamics using Euler integration
    qdAct = qdAct + dt * (M \ (Torque - C)')'; % Update velocity
    qAct = qAct + dt * qdAct; % Update position
    
    % Store for visualization
    J_current = jacobian_2link(qAct(1), qAct(2), l1, l2);
    F = (-K * (inv(J_current') * Torque'))';
    tau(i, :) = Torque'; % Save torque
    Force(i, :) = F; % Save force
    xPos(i, :) = forward_kinematics(qAct(1), qAct(2), l1, l2); % Save position
    qPos(i,:) = qAct; % Save joint angles
    
    % Update qd_Des_Prev for the next iteration
    qd_Des_Prev = qdAct;
end


% Plot trajectory
figure(1); hold on; grid on;
plot(t, qDes(:,1), '--', 'DisplayName', 'q1 Desired');
plot(t, qDes(:,2), '--', 'DisplayName', 'q2 Desired');
plot(t, qPos(:,1), '-', 'DisplayName', 'q1 Actual');
plot(t, qPos(:,2), '-', 'DisplayName', 'q2 Actual');
xlabel('Time (s)');
ylabel('Joint Angles (rad)');
title('Joint Space Trajectory Following');
legend show;

figure(2); hold on; grid on;
plot(qDes(1, 1), qDes(2, 1), '*', 'DisplayName', 'Equilibrium Point');
plot(qPos(1, 1), qPos(2, 1), '*', 'DisplayName', 'Starting Point');
quiver(qPos(:,1), qPos(:,2), tau(:,1), tau(:, 2));
legend show;

% Cartesian Space
figure(3); hold on; grid on;
plot(t, xPos(:,1), '--', 'DisplayName', 'EndEffector Desired X');
plot(t, xPos(:,2), '--', 'DisplayName', 'EndEffector Desired Y');
xlabel('Time (s)');
ylabel('End Effector ()');
title('Cartesian Space Trajectory Following');
legend show;


figure(4); hold on; grid on;
step = 50:501;
plot(xPos(1,1),xPos(1,2),'*', 'DisplayName', 'First Position')
plot(xPos(end,1),xPos(end,2),'*', 'DisplayName', 'Last Position')
quiver(xPos(step,1), xPos(step,2), Force(step, 1), Force(step, 2));
legend show;






