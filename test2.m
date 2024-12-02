close all;

% Parameters
l1 = 1; l2 = 1;          % Link lengths
m1 = 1; m2 = 1;          % Masses
g = -9.81;               % Gravity
F_end_effector = [0; 0]; % External force [Fx; Fy]
K = 1;                   % Stiffness gain

% Desired joint position 
xDes = 0.75; yDes = 0.75;
[q1, q2] = inverse_kinematics(xDes, yDes, l1, l2);
qDes = [q1 q2];  % Desired joint angles

% Actual Cartesian positions (meshgrid)
xAct = 1.2; yAct = 1.2;
[q1, q2] = inverse_kinematics(xAct, yAct, l1, l2);
qAct = [q1 q2];

% Calculate the error between the desired and actual joint positions
error = qDes - qAct;  % Column vector for error

% Calculate the Jacobian at the current joint angles
J_current = jacobian_2link(q1, q2, l1, l2);

% Compute the force using the pseudo-inverse of the Jacobian
Force = (-K * (inv(J_current') * error'))';  % Force vector (row)


%%%%%%%%
% Time
% Parameters
T = 10; % Total time duration (seconds)
x0 = T / 2; % Midpoint of the transition
k = 5; % Steepness of the transition

% Time vector (from 0 to T)
t = linspace(0, T, 1000);

% Sigmoid function
time_step = 1 ./ (1 + exp(-k * (t - x0)));

% Derivative of the sigmoid function
time_impulse = (k * exp(-k * (t - x0))) ./ (1 + exp(-k * (t - x0))).^2;

% Step and Impulse Force Fields
Step_ForceField = Force' * time_step;
Impulse_ForceField = Force' * time_impulse;

% Plot joint angles (qDes vs qAct)
figure;

% Plot the workspace with desired and actual positions
subplot(2, 2, 1);
hold on;
plot([0, l1*cos(qDes(1)) + l2*cos(qDes(1) + qDes(2))], [0, l1*sin(qDes(1)) + l2*sin(qDes(1) + qDes(2))], 'r', 'LineWidth', 2);
plot([0, l1*cos(qAct(1)) + l2*cos(qAct(1) + qAct(2))], [0, l1*sin(qAct(1)) + l2*sin(qAct(1) + qAct(2))], 'b', 'LineWidth', 2);
title('Arm Configurations in Workspace');
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend('Desired Position', 'Actual Position');
axis equal;
grid on;

% Plot the Force field for Step_ForceField
subplot(2, 2, 2);
quiver(0, 0, Step_ForceField(1), Step_ForceField(2), 0, 'LineWidth', 2, 'MaxHeadSize', 1);
axis equal;
xlim([-2, 2]);
ylim([-2, 2]);
title('Step Force Field');
xlabel('Force in X (N)');
ylabel('Force in Y (N)');
grid on;

% Plot the Force field for Impulse_ForceField
subplot(2, 2, 3);
quiver(0, 0, Impulse_ForceField(1), Impulse_ForceField(2), 0, 'LineWidth', 2, 'MaxHeadSize', 1);
axis equal;
xlim([-2, 2]);
ylim([-2, 2]);
title('Impulse Force Field');
xlabel('Force in X (N)');
ylabel('Force in Y (N)');
grid on;

% Plot the time series of the step and impulse force fields
subplot(2, 2, 4);
plot(t, Step_ForceField(1, :), 'r', 'LineWidth', 2);
hold on;
plot(t, Step_ForceField(2, :), 'b', 'LineWidth', 2);
plot(t, Impulse_ForceField(1, :), 'g', 'LineWidth', 2);
plot(t, Impulse_ForceField(2, :), 'k', 'LineWidth', 2);
title('Step and Impulse Force Fields Over Time');
xlabel('Time (s)');
ylabel('Force (N)');
legend('Step Force X', 'Step Force Y', 'Impulse Force X', 'Impulse Force Y');
grid on;
