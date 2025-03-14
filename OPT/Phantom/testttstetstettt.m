clear; clc;

% Define desired trajectory and Middle Points
qDes = [0.1914, -0.0445, 0.3336];

qMid = zeros(5,3);
qMid(1,:) = IK(0.025, 0, 0.01);
qMid(2,:) = IK(0.03,  0, 0.015);
qMid(3,:) = IK(0.035, 0, 0.02);
qMid(4,:) = IK(0.04,  0, 0.025);
qMid(5,:) = IK(0.045, 0, 0.03);

% Parameters
time = 10;  % Time
wn = [1 1 10];  % Prefilter Omega     
kj = [20 20 20];  % Spring constants
bj = [5 5 5];  % Damping constants
wt = [400, 0.001, 10000];  % Weights [qDes, Time, qMid]

% Optimization setup

[initT, initY] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj), [0 time], zeros(12, 1));
[xInit, yInit, zInit] = FK(initY(:,7), initY(:,8), initY(:,9));
[xD, yD, zD] = FK(initY(:,1), initY(:,2), initY(:,3));

% Prepare the figure for animation

qDes = 0.1914 * ones(size(initT));
% figure(1); hold on; grid on;
% % plot([0 initT'],[0 qDes(:,1)'])
% plot(initT,initY(:,1))
% text(4, 0.3, 'Under-Damped System', 'FontSize', 12, 'Color', 'r')
% text(2,.05, 'Over-damped System', 'FontSize', 12, 'Color', 'r')
% 
% text(3,.13, 'Critical damping System', 'FontSize', 12, 'Color', 'r')
% 
% 
% xlabel('Time(s)')
% ylabel('Joint Position(rad)')
% title('Impact of Sigma on Second-Order Filter Dynamics')
% 

% figure(1); hold on; grid on;



figure(2); hold on; grid on;
subplot(3,1,1)
plot(initT, initY(:,7), 'LineWidth', 2)
xlabel('Time (s)', 'FontSize', 16, 'FontName', 'Times New Roman')
ylabel('Joint 1 (rad)', 'FontSize', 16, 'FontName', 'Times New Roman')
title('First Joint Position', 'FontSize', 20, 'FontName', 'Courier New') 
grid on

subplot(3,1,2)
plot(initT, initY(:,9), 'LineWidth', 2)
xlabel('Time (s)', 'FontSize', 16, 'FontName', 'Times New Roman')
ylabel('Joint 2 (rad)', 'FontSize', 16, 'FontName', 'Times New Roman')
title('Second Joint Position', 'FontSize', 20, 'FontName', 'Courier New') 
grid on;

subplot(3,1,3)
plot(xInit, zInit, '-', 'LineWidth', 2)
grid on
xlabel('X Axis', 'FontSize', 16, 'FontName', 'Times New Roman')
ylabel('Y Axis', 'FontSize', 16, 'FontName', 'Times New Roman')
title('Cartesian Space Trajectory', 'FontSize', 20, 'FontName', 'Courier New') 

% myTwolinkwithprefilter function
function dxdt = myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj)
    zeta = 3;
    A = [zeros(3,3) eye(3);
        -eye(3)*diag(wn).^2  -eye(3)*2*zeta*diag(wn)];
    B = [zeros(3,3); diag(wn).^2];
  
    q = x(7:9);
    qd = x(10:12);

    Kp = -diag(kj);  
    Kd = -diag(bj);  

    [M, C, G] = compute_M_C_G(q(1), q(2), q(3), qd(1), qd(2), qd(3));
    torque = Kd * (qd - x(4:6)) + Kp * (q - x(1:3));
    qdd = M \ (C * qd + torque);

    dxdt = [A * x(1:6) + B * qDes(:); qd; qdd];
end

function [x, y, z] = FK(q1, q2, q3)
    l1 = 0.208; 
    l2 = 0.168;  
    x = sin(q1) .* (l1 * cos(q2) + l2 * sin(q3));
    y = l2 - l2 * cos(q3) + l1 * sin(q2);
    z = -l1 + cos(q1) .* (l1 * cos(q2) + l2 * sin(q3));
end

function [q1, q2, q3] = IK(x, y, z)
    l1 = 0.208; 
    l2 = 0.168;  
    q1 = atan2(x, z + l1);

    R = sqrt(x^2 + (z + l1)^2);
    r = sqrt(x^2 + (y - l2)^2 + (z + l1)^2);
    Beta  = atan2(y - l2, R);
    Gamma = acos((l1^2 + r^2 - l2^2) / (2 * l1 * r));
    q2 = Gamma + Beta;

    Alpha = acos((l1^2 + l2^2 - r^2) / (2 * l1 * l2));
    q3 = q2 + Alpha - pi/2;
end

% clear; clc;
% 
% % Define desired trajectory and Middle Points
% qDes = [0.1914, -0.0445, 0.3336];
% 
% qMid = zeros(5,3);
% qMid(1,:) = IK(0.025, 0, 0.01);
% qMid(2,:) = IK(0.03,  0, 0.015);
% qMid(3,:) = IK(0.035, 0, 0.02);
% qMid(4,:) = IK(0.04,  0, 0.025);
% qMid(5,:) = IK(0.045, 0, 0.03);
% 
% % Parameters
% time = 10;  % Time
% wn = [1 1 1];  % Prefilter Omega     
% kj = [20 20 20];  % Spring constants
% bj = [5 5 5];  % Damping constants
% wt = [400, 0.001, 10000];  % Weights [qDes, Time, qMid]
% 
% % Optimization setup
% initParams = [time wn bj kj]; % Initial guess
% 
% [init_T, init_Y] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj), [0 time], zeros(12, 1));
% [xInit, yInit, zInit] = FK(init_Y(:,7), init_Y(:,8), init_Y(:,9));
% [xD, yD, zD] = FK(init_Y(:,1), init_Y(:,2), init_Y(:,3));
% 
% figure(1); hold on; grid on;
% 
% subplot(3,1,1)
% plot(init_T,init_Y(:,7), 'LineWidth', 2)
% xlabel('Time (s)', 'FontSize', 16, 'FontName', 'Times New Roman')
% ylabel('Joint 1 (rad)', 'FontSize', 16, 'FontName', 'Times New Roman')
% title('First Joint Position', 'FontSize', 20, 'FontName', 'Courier New') % Changed font
% grid on
% 
% subplot(3,1,2)
% plot(init_T,init_Y(:,9), 'LineWidth', 2)
% xlabel('Time (s)', 'FontSize', 16, 'FontName', 'Times New Roman')
% ylabel('Joint 2 (rad)', 'FontSize', 16, 'FontName', 'Times New Roman')
% title('Second Joint Position', 'FontSize', 20, 'FontName', 'Courier New') % Changed font
% grid on;
% 
% subplot(3,1,3) 
% plot(xInit, zInit, '-', 'LineWidth', 2)
% grid on
% xlabel('X Axis', 'FontSize', 16, 'FontName', 'Times New Roman')
% ylabel('Y Axis', 'FontSize', 16, 'FontName', 'Times New Roman')
% title('Cartesian Space Trajectory', 'FontSize', 20, 'FontName', 'Courier New') % Changed font
% 
% 
% 
% % myTwolinkwithprefilter function
% function dxdt= myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj)
%     zeta = 1;
%     A = [zeros(3,3) eye(3);
%         -eye(3)*diag(wn).^2  -eye(3)*2*zeta*diag(wn)];
%     B = [zeros(3,3); diag(wn).^2];
% 
%     q   = x(7:9);
%     qd  = x(10:12);
% 
%     Kp = -diag(kj);  
%     Kd = -diag(bj);  
% 
%     [M, C, G] = compute_M_C_G(q(1), q(2), q(3), qd(1), qd(2), qd(3));
%     torque = Kd * (qd - x(4:6)) + Kp * (q - x(1:3));
%     qdd = M \ (C * qd + torque);
% 
%     dxdt = [A*x(1:6) + B*qDes(:); qd; qdd];
% 
% end
% 
% function [x, y, z] = FK(q1, q2, q3)
%     l1 = 0.208; 
%     l2 = 0.168;  
%     x = sin(q1) .* (l1 * cos(q2) + l2 * sin(q3));
%     y = l2 - l2 * cos(q3) + l1 * sin(q2);
%     z = -l1 + cos(q1) .* (l1 * cos(q2) + l2 * sin(q3));
% end
% 
% function [q1, q2, q3] = IK(x, y, z)
%     l1 = 0.208; 
%     l2 = 0.168;  
%     q1 = atan2(x, z + l1);
% 
%     R = sqrt(x^2 + (z + l1)^2);
%     r = sqrt(x^2 + (y - l2)^2 + (z + l1)^2);
%     Beta  = atan2(y - l2, R);
%     Gamma = acos((l1^2 + r^2 - l2^2) / (2 * l1 * r));
%     q2 = Gamma + Beta;
% 
%     Alpha = acos((l1^2 + l2^2 - r^2) / (2 * l1 * l2));
%     q3 = q2 + Alpha - pi/2;
% end
