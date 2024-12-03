close all
clear

% Parameters
l1 = 1; l2 = 1; % Link lengths
m1 = 1; m2 = 1; % Masses
g = 0;

qDes = [pi/8; 0];
qdDes = [0;0];

% Define LQR parameters
Q = diag([20, 20, 20, 20]); % State cost matrix (adjust to tune)
R = diag([.01, .01]);         % Control input cost matrix (adjust to tune)

% Solve system dynamics
odefun = @(t,x) mysf_with_lqr(t,x,l1,l2,m1,m2,qDes,qdDes,Q,R);
[t, y] = ode45(odefun, [0 20], [0; 0; 0; 0]);

% Create a linearly spaced time vector
t_uniform = linspace(0, 20, 20000); % Adjust resolution as needed

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
    
    % Compute end-effector position
    xPos(i, :) = forward_kinematics(y_interp(i, 1), y_interp(i, 2), l1, l2);
end

% % Plot the results
figure(1); hold on; grid on;
plot(t_uniform, y_interp(:,1)); title('Joint 1 Position');
xlabel('Time (s)'); ylabel('Position (rad)');

figure(2); hold on; grid on;
plot(t_uniform, y_interp(:,2)); title('Joint 2 Position');
xlabel('Time (s)'); ylabel('Position (rad)');

% --- Function Definitions ---
function dx = mysf_with_lqr(t, x, l1, l2, m1, m2, qDes, qdDes, Q, R)
    % Extract state variables
    qAct = x(1:2);      % Joint positions
    qdAct = x(3:4);     % Joint velocities
    
    % Compute dynamics
    M = mass_matrix(qAct(1), qAct(2), l1, l2, m1, m2);
    G = gravity_vector(qAct(1), qAct(2), l1, l2, m1, m2, 0, 0);
    C = coriolis_matrix(qAct(1), qAct(2), qdAct(1), qdAct(2), l1, l2, m1, m2);

    % State error
    e = [qAct - qDes; qdAct - qdDes];
    
    % Linearized dynamics (A and B matrices)
    A = [zeros(2), eye(2);
         eye(2), zeros(2)];
    B = [zeros(2); eye(2)];
    
    % LQR gain matrix
    [K,S,P] = lqr(A, B, Q, R);
    
    % Control input using LQR
    Torque = -K * e;
    
    % Compute accelerations
    qddAct = M \ (Torque - C(:));
    
    % Return derivatives of state variables
    dx = [qdAct(:); qddAct(:)];
end

