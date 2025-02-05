% Parameters
l1 = 1; l2 = 1; % Link lengths
m1 = 1; m2 = 1; % Masses
g = 0;

tspan = [0 20]; % Total simulation time

% Prefilter (Low-Pass)
Wn = 2;
sigma = 1;
num = Wn^2;
den = [1 2*sigma*Wn Wn^2];
Prefilter = tf(num,den);

% Desired Joint Trajectory
q = [    -0.4240,   2.4189; % MAIN

    % -0.2936    2.3392 ;% second middle
    % -0.2       2.26;
    % -0.1205    2.2065 ;% first middle
    % -0.04      2.13;
    % 0.0337     2.0601; % second middle
    % 0.07       2;
    0.1296    1.9552; % MAIN
    % 0.1       1.92;
    % 0.0821     1.8965 ;% second middle
    % 0.055      1.84;
    % 0.0316     1.7913 ;% first middle
    % 0.015      1.72;
    % 0.0050     1.6659 ;% second middle
    % 0.002      1.62;

    0.0,       1.5708; % MAIN
    % -0.1002    1.6659;% second middle
    % -0.2522    1.7913 ;% first middle
    % -0.4078    1.8965;% second middle
    -0.5139   1.9552; % MAIN
    % -0.5229    2.0601 ;%  second middle
    % -0.5153    2.2065 ;% first middle
    % -0.4749    2.3392; %  second middle
    -0.4240,   2.4189 % MAIN
    ]; % Modify as needed

% qDes, qDesF and change the first 300 value with the home  position
qDes = repelem(q, 400, 1);
tNorm = linspace(0, tspan(2), length(qDes)); % Time vector
tFilt = linspace(0, tspan(2), length(qDes)); % Time vector

qDesF = [lsim(Prefilter, qDes(:,1), tNorm, qDes(1,1)), lsim(Prefilter, qDes(:,2), tNorm, qDes(1,2))];
step = 300;
qDesF(1:step,:) = [ones(step,1) * -0.4240, ones(step,1)*2.4189];  % starts from [0,0]

xDesNorm = forward_kinematics(qDes(:, 1), qDes(:, 2), l1, l2);
xDesFilt = forward_kinematics(qDesF(:, 1), qDesF(:, 2), l1, l2);

% Controller gains
K = 80; % Proportional gain
B = 60; % Derivative gain

% Initial joint angles and velocities
q0 = [-0.4240; 2.4189]; % Initial joint angles
qd0 = [0; 0]; % Initial joint velocities
initial_state = [q0; qd0];

% Solve the ODE
[tNormOut, stateNorm] = ode45(@(tNorm, x) robot_dynamics(tNorm, x, l1, l2, m1, m2, g, K, B, qDes, tspan), tspan, initial_state);
[tFiltOut, stateFilt] = ode45(@(tFilt, x) robot_dynamics(tFilt, x, l1, l2, m1, m2, g, K, B, qDesF, tspan), tspan, initial_state);

% Extract results
qPosNorm = stateNorm(:, 1:2); % Joint positions
qVelNorm = stateNorm(:, 3:4); % Joint velocities
qPosFilt = stateFilt(:, 1:2); % Joint positions
qVelFilt = stateFilt(:, 3:4); % Joint velocities

xPosActNorm = forward_kinematics(qPosNorm(:, 1), qPosNorm(:, 2), l1, l2);
xPosActFilt = forward_kinematics(qPosFilt(:, 1), qPosFilt(:, 2), l1, l2);

xVelActJNorm = zeros(size(xPosActNorm));
xVelActJFilt = zeros(size(xPosActFilt));

for i=1:length(qPosNorm)
    xVelActJNorm(i,:) = jacobian_2link(qPosNorm(i,1),qPosNorm(i,2),1,1) * qVelNorm(i,:)';
end
for i=1:length(qPosFilt)
    xVelActJFilt(i,:) = jacobian_2link(qPosFilt(i,1),qPosFilt(i,2),1,1) * qVelFilt(i,:)';
end


figure(4); hold on; grid on;
plot(xDesNorm(:, 1), xDesNorm(:, 2), '-', 'DisplayName', 'End Effector Desired - NORMAL');
plot(xPosActFilt(:, 1), xPosActFilt(:, 2), '-', 'DisplayName', 'End Effector Actual Position - Filtered');
xlabel('X-Axis')
ylabel('Y-Axis')
legend show
%%

% Plot trajectory
figure(1); hold on; grid on;
plot(tNorm, qDes(:,1))
plot(tFilt, qDesF(:,1))
xlabel('Time (s)')
ylabel('Joint Position (rad)')
legend('Desired Joint1 Position','Filtered Joint1 Position')

% figure(2); hold on; grid on;
% plot(tNorm, qDes(:,2))
% plot(tFilt, qDesF(:,2))
% xlabel('Time (s)')
% ylabel('Joint Position (rad)')
% legend('Desired Joint2 Position','Filtered Joint2 Position')

figure(3); hold on; grid on;
plot(xDesNorm(:, 1), xDesNorm(:, 2), '-', 'DisplayName', 'End Effector Desired - NORMAL');
plot(xDesFilt(:, 1), xDesFilt(:, 2), '-', 'DisplayName', 'End Effector Desired - FILTERED');
xlabel('X-Axis')
ylabel('Y-Axis')
legend show

figure(4); hold on; grid on;
plot(xDesNorm(:, 1), xDesNorm(:, 2), '-', 'DisplayName', 'End Effector Desired - NORMAL');
plot(xPosActFilt(:, 1), xPosActFilt(:, 2), '-', 'DisplayName', 'End Effector Actual Position - Filtered');
xlabel('X-Axis')
ylabel('Y-Axis')
legend show

% 
% figure(4); hold on; grid on;
% plot(xPosActNorm(:, 1), xPosActNorm(:, 2), '-', 'DisplayName', 'End Effector Actual Position - Normal');
% plot(xPosActFilt(:, 1), xPosActFilt(:, 2), '-', 'DisplayName', 'End Effector Actual Position - Filtered');
% xlabel('X-Axis')
% ylabel('Y-Axis')
% legend show

figure(5); hold on; grid on;
plot(tNormOut,xVelActJNorm(:, 1:2), '-', 'DisplayName', 'End Effector Actual Velocity 1 - Normal');
plot(tFiltOut,xVelActJFilt(:, 1:2), '-', 'DisplayName', 'End Effector Actual Velocity 1 - Filtered');
xlabel('Time')
ylabel('vel')
legend show



figure(6); hold on; grid on;
plot(tNormOut, qPosNorm(:,1))
plot(tFiltOut, qPosFilt(:,1))
legend('Normal joint1 position','Filtered joint1 position')

% figure(7); hold on; grid on;
% plot(tNormOut, qPosNorm(:,2))
% plot(tFiltOut, qPosFilt(:,2))
% legend('Normal joint2 position','Filtered joint2 position')

figure(8); hold on; grid on;
plot(tNormOut, qVelNorm(:,1))
plot(tFiltOut, qVelFilt(:,1))
legend('Normal joint1 velocity','Filtered joint1 velocity')

% figure(9); hold on; grid on;
% plot(tNormOut, qVelNorm(:,2))
% plot(tFiltOut, qVelFilt(:,2))
% legend('Normal joint2 velocity','Filtered joint2 velocity')


figure(10); hold on; grid on;
plot(tNormOut,xPosActNorm(:,1))
plot(tFiltOut,xPosActFilt(:,1))
legend('Output position Normal','Output position filtered')

figure(11); hold on; grid on;
plot(tNormOut,xVelActJNorm(:,1))
plot(tFiltOut,xVelActJFilt(:,1))
legend('Output velocity Normal','Output velocity filtered')


% Functions
function dxdt = robot_dynamics(t, x, l1, l2, m1, m2, g, K, B, Q, tspan)
% Unpack state variables
q = x(1:2);
qd = x(3:4);


% Determine current waypoint
num_waypoints = size(Q, 1);
val = tspan(2) / num_waypoints;
waypoint_index = min(floor(t / val) + 1, num_waypoints);
q_des = Q(waypoint_index, :)';

% Desired velocity
qd_des = [0; 0];

% Errors
e = q - q_des;
eDot = qd - qd_des;

% Dynamics
M = mass_matrix(q(1), q(2), l1, l2, m1, m2);
G = gravity_vector(q(1), q(2), l1, l2, m1, m2, g, g);
C = coriolis_matrix(q(1), q(2), qd(1), qd(2), l1, l2, m1, m2);

% PD control with gravity compensation
Torque = -K * e - B * eDot;

% State derivatives
qdd = M \ (Torque - C * qd - G);
dxdt = [qd; qdd];
end


function P = forward_kinematics(q1, q2, l1, l2)
x = l1 * cos(q1) + l2 * cos(q1 + q2);
y = l1 * sin(q1) + l2 * sin(q1 + q2);
P = [x, y];
end

function P = FK(q1, q2, l1, l2)
L = [ 0 0 ; l1 0; l1 l2];
C = [cos(q1) sin(q1);
     cos(q1+q2) sin(q1+q2)];
P = L * C;

end

function qDes = inverse_kinematics(x, y, l1, l2)
r = sqrt(x^2 + y^2);
if r > (l1 + l2) || r < abs(l1 - l2)
    error('Target point is out of reach');
end
cos_q2 = (r^2 - l1^2 - l2^2) / (2 * l1 * l2);
sin_q2 = sqrt(1 - cos_q2^2);
q2 = atan2(sin_q2, cos_q2);
phi = atan2(y, x);
psi = atan2(l2 * sin(q2), l1 + l2 * cos(q2));
q1 = phi - psi;
qDes = [q1, q2];
end


function M = mass_matrix(q1, q2, l1, l2, m1, m2)
M = [
    (m1 + m2) * l1^2 + m2 * l2^2 + 2 * m2 * l1 * l2 * cos(q2), m2 * l2^2 + m2 * l1 * l2 * cos(q2);
    m2 * l2^2 + m2 * l1 * l2 * cos(q2), m2 * l2^2
    ];
end


function G = gravity_vector(q1, q2, l1, l2, m1, m2, g1, g2)
G = [
    -(m1 + m2) * g1 * l1 * sin(q1) - m2 * g2 * l2 * sin(q1 + q2);
    -m2 * g2 * l2 * sin(q1 + q2)
    ];
end


function C = coriolis_matrix(q1, q2, q1d, q2d, l1, l2, m1, m2)
h = -m2 * l1 * l2 * sin(q2);
C = [
    h * q2d, h * (q1d + q2d);
    -h * q1d, 0
    ];
end
