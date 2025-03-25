clc; clear; close all;

% Define system parameters for each joint
wn = [5, 8, 10];         % Natural frequencies for joints (rad/s)
zeta = [0.7, 1.0, 0.5];  % Damping ratios for joints

% Number of joints
n_joints = 3;

% Create state-space models for each joint
sys = cell(1, n_joints);
for i = 1:n_joints
    A = [0 1 0; 0 0 1; -wn(i)^3 -3*zeta(i)*wn(i)^2 -3*zeta(i)^2*wn(i)];
    B = [0; 0; wn(i)^3];
    C = [1 0 0]; % Output is the joint angle θ
    D = 0;
    
    % Create state-space system
    sys{i} = ss(A, B, C, D);
end

% Step Response for Each Joint
figure;
hold on;
for i = 1:n_joints
    step(sys{i});
end
title('Step Response for Each Joint');
legend('Joint 1', 'Joint 2', 'Joint 3');
xlabel('Time (s)');
ylabel('Response');
grid on;
hold off;

% Bode Plots for Each Joint
figure;
hold on;
for i = 1:n_joints
    bode(sys{i});
end
title('Bode Plot for Each Joint');
legend('Joint 1', 'Joint 2', 'Joint 3');
grid on;
hold off;

