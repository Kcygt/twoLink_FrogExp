close all;

% Parameters
l1 = 1; l2 = 1;          % Link lengths
m1 = 1; m2 = 1;          % Masses
g = -9.81;               % Gravity
F_end_effector = [0; 0]; % External force [Fx; Fy]
K = 1;                 % Stiffness gain

% Desired joint position 
xDes = 0.75; yDes = 0.75;
qDes = inverse_kinematics(xDes, yDes, l1, l2);
 
% Desired Cartesian positions (meshgrid)
xAct = 1.2; yAct = 1.2;
qAct = inverse_kinematics(xAct, yAct, l1, l2);

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
k = 1; % Steepness of the transition

% Time vector (from 0 to T)
t = linspace(0, T, 10);

% Sigmoid function
time_step= 1 ./ (1 + exp(-k * (t - x0)));

% Derivative of the sigmoid function
time_impulse = (k * exp(-k * (t - x0))) ./ (1 + exp(-k * (t - x0))).^2;

% Calculted Forces
Step_ForceField = Force' * time_step;
Impulse_ForceField = Force' * time_impulse;


figure(1); hold on; grid on;
plot(qDes(1),qDes(2),'*',qAct(1),qAct(2),'o')
legend('Desired Joint Position','Actual Joint Position')

