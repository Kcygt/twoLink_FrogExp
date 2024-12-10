close all;

% Parameters
l1 = 1; l2 = 1;           % Link lengths
m1 = 1; m2 = 1;           % Masses
g = 0;                    % Gravity
F_end_effector = [0; 0];  % External force [Fx; Fy]
K = 100;                  % Stiffness gain
qAct = [];                % Actual joint positions (initially empty)

% Desired joint position 
xDes = [0.75; 0.75];
[q1, q2] = inverse_kinematics(xDes(1), xDes(2), l1, l2);
qDes = [q1 q2];  % Desired joint angles

% Desired Cartesian positions (meshgrid)
[xAct, yAct] = meshgrid(0.3:0.3:1.2, 0.3:0.3:1.2);
xPos = [xAct(:), yAct(:)];

% Preallocate force array
Force = zeros(numel(xAct), 2);  % Size based on meshgrid
tau = zeros(numel(xAct), 2);  % Size based on meshgrid
TAU = zeros(numel(xAct), 2);  % Size based on meshgrid

% Loop through the meshgrid to calculate actual joint positions and forces
for ii = 1:numel(xAct)
    [q1, q2] = inverse_kinematicsE(xAct(ii), yAct(ii), l1, l2);
    qAct = [qAct; q1, q2];  % Store joint angles as rows
    
    % Calculate the error between the desired and actual joint positions
    e = qDes' - [q1; q2];  % Column vector for error
    
    % Calculate the Jacobian at the current joint angles
    J_current = jacobian_2link(q1, q2, l1, l2);

    % Compute the force using the pseudo-inverse of the Jacobian
    Force(ii, :) = (K * (inv(J_current') * e))';  % Force vector (row)
    tau(ii,:) = K * e ;
    TAU(ii, :) = K * e * exp(-e'*K*e);
end

% Plot joint and end-effector positions
figure(1); hold on; grid on;
plot(qDes(1), qDes(2), '*', 'DisplayName', 'Desired');
plot(qAct(:, 1), qAct(:, 2), 'o', 'DisplayName', 'Actual');
quiver(qAct(:, 1), qAct(:, 2), tau(:, 1), tau(:, 2), 'DisplayName', 'Force');
xlabel('Joint Angle 1 (rad)');
ylabel('Joint Angle 2 (rad)');
title('Actual and Desired Joint Positions with Torques');
legend show;

% Plot joint and end-effector positions
figure(2); hold on; grid on;
plot(xDes(1), xDes(2), '*', 'DisplayName', 'Desired');
plot(xPos(:, 1), xPos(:, 2), 'o', 'DisplayName', 'Actual');
quiver(xPos(:, 1), xPos(:, 2), Force(:, 1), Force(:, 2), 'DisplayName', 'Force');
xlabel('Joint Angle 1 (rad)');
ylabel('Joint Angle 2 (rad)');
title('Actual and Desired Cartesian Positions with Forces');
legend show;


% Function Definitions
function [x, y] = forward_kinematics(q1, q2, l1, l2)
    % Calculate end-effector position from joint angles
    x = l1*cos(q1) + l2*cos(q1 + q2);
    y = l1*sin(q1) + l2*sin(q1 + q2);
end

function [q1, q2] = inverse_kinematics(x, y, l1, l2)
    % Inverse kinematics for a 2-link manipulator
    r = sqrt(x^2 + y^2);  % Distance to target
    if r > (l1 + l2) || r < abs(l1 - l2)
        error('Target point is out of reach');
    end
    cos_q2 = (r^2 - l1^2 - l2^2) / (2 * l1 * l2);
    sin_q2 = sqrt(1 - cos_q2^2);
    q2 = atan2(sin_q2, cos_q2);
    phi = atan2(y, x);
    psi = atan2(l2 * sin(q2), l1 + l2 * cos(q2));
    q1 = phi - psi;
end

function J = jacobian_2link(q1, q2, l1, l2)
    % Compute the Jacobian matrix for a 2-link manipulator
    J = [
        -l1*sin(q1) - l2*sin(q1 + q2), -l2*sin(q1 + q2);
         l1*cos(q1) + l2*cos(q1 + q2),  l2*cos(q1 + q2)
    ];
end

function [q1, q2] = inverse_kinematicsE(x, y, l1, l2)
    % Inverse kinematics with error handling for unreachable targets
    r = sqrt(x^2 + y^2);
    if r > (l1 + l2) || r < abs(l1 - l2)
        warning('Target point (%f, %f) is out of reach', x, y);
        q1 = NaN; q2 = NaN;  % Return NaN for unreachable points
        return;
    end
    cos_q2 = (r^2 - l1^2 - l2^2) / (2 * l1 * l2);
    sin_q2 = sqrt(1 - cos_q2^2); 
    q2 = atan2(sin_q2, cos_q2);
    phi = atan2(y, x);
    psi = atan2(l2 * sin(q2), l1 + l2 * cos(q2));
    q1 = phi - psi;
end
