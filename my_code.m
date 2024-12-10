close all
clear

% Parameters
l1 = 1; l2 = 1; % Link lengths
m1 = 1; m2 = 1; % Masses
g = 0;

qDes = [pi/8; pi/4];
qdDes = [0;0];

% Controller gains
K = 2; % Proportional gain (adjust as needed)
B = 5; % Derivative gain (adjust as needed)

% Solve system dynamics
odefun = @(t,x) mysf(t,x,l1,l2,m1,m2,K,B,qDes,qdDes);
[t, y] = ode45(odefun, [0 20], [0; 0; 0; 0]);

% Create a linearly spaced time vector
t_uniform = linspace(0, 20, 10000); % Adjust resolution as needed

% Interpolate the results
y_interp = interp1(t, y, t_uniform);

% Preallocate arrays for forces, positions, and torques
F = zeros(length(t_uniform), 2);
xPos = zeros(length(t_uniform), 2);
Torque = zeros(length(t_uniform), 2);

% Process interpolated results
for i = 1:length(t_uniform)

    % Jacobian matrix at the current configuration
    J_current = jacobian_2link(y_interp(i, 1), y_interp(i, 2), l1, l2);
    
    % Compute position and velocity errors
    e = [y_interp(i, 1); y_interp(i, 2)] - [pi/4; 0];
    eDot = [y_interp(i, 3); y_interp(i, 4)] - [0; 0];
    
    % PD control law for torque
    % Torque(i, :) = (K * -e + B * -eDot)';
    Torque(i, :) = K * -e * exp(-e'*K*e) + B * -eDot;

    % Convert torque to Cartesian forces using the Jacobian
    F(i, :) = (-K * (inv(J_current') * Torque(i, :)'))';

    % Compute end-effector position
    xPos(i, :) = forward_kinematics(y_interp(i, 1), y_interp(i, 2), l1, l2);
end

% % Plot the results
figure(1); hold on; grid on;
quiver(y_interp(:,1),y_interp(:,2),Torque(:,1),Torque(:,2))

figure(2); hold on; grid on;
quiver(xPos(:,1),xPos(:,2),F(:,1),F(:,2))

figure(3); hold on; grid on;
plot(y_interp(:,1))

% --- Function Definitions ---

function dx = mysf(t, x, l1, l2, m1, m2, K, B,qDes,qdDes)
    % Extract state variables
    qAct = x(1:2);      % Joint positions
    qdAct = x(3:4);     % Joint velocities
    % Compute dynamics
    M = mass_matrix(qAct(1), qAct(2), l1, l2, m1, m2);
    G = gravity_vector(qAct(1), qAct(2), l1, l2, m1, m2, 0, 0);
    C = coriolis_matrix(qAct(1), qAct(2), qdAct(1), qdAct(2), l1, l2, m1, m2);
    
    % PD control law
    e = qAct - qDes;
    eDot = qdAct - qdDes;
    Torque = K * -e * exp(-e'*K*e) + B * -eDot;
    % Torque = K * -e + B * -eDot;

    % Compute accelerations
    qddAct = M \ (Torque - C(:));
    
    % Return derivatives of state variables
    dx = [qdAct(:); qddAct(:)];
end
