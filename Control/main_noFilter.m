% Parameters
l1 = 1; l2 = 1; % Link lengths
m1 = 1; m2 = 1; % Masses
g = 0;

tspan = [0 20]; % Total simulation time

% Desired Joint Trajectory
qDes = [ -0.4240,   2.4189;  
          0.1296    1.9552;  
          0.0,      1.5708;  
         -0.5139    1.9552;  
         -0.4240,   2.4189  
        ]; 

% Controller gains
K = 80; % Proportional gain
B = 40; % Derivative gain

% Initial joint angles and velocities
q0 = [-0.4240; 2.4189]; % Initial joint angles
qd0 = [0; 0]; % Initial joint velocities
initial_state = [q0; qd0];

[tOut, state] = ode45(@(t, x) robot_dynamics(t, x, l1, l2, m1, m2, g, K, B, qDes), tspan, initial_state);
qAct = state(:, 1:2); % Joint positions
qdAct = state(:, 3:4); % Joint velocities

xAct= forward_kinematics(qAct(:, 1), qAct(:, 2), l1, l2);

xdAct= zeros(size(xAct));
for i=1:length(qAct)
    xdAct(i,:) = jacobian_2link(qAct(i,1),qAct(i,2),1,1) * qdAct(i,:)';
end


% Visualize Robot Movement
figure(1); hold on; axis equal; grid on;
xlabel('X-axis'); ylabel('Y-axis');
title('2-Link Robot Movement');
xlim([-2, 2]); ylim([-2, 2]);

% Number of timesteps
num_timesteps = size(qAct, 1);

% Plot setup
robot_plot = plot(0, 0, 'o-', 'LineWidth', 2, 'MarkerSize', 6);

% Animate the robot
for t = 1:num_timesteps
    % Get joint angles at timestep t
    q1 = qAct(t, 1);
    q2 = qAct(t, 2);
    
    % Compute forward kinematics
    P = FK(q1, q2, l1, l2);
    
    % Update robot plot
    set(robot_plot, 'XData', [0, P(2, 1), P(3, 1)], ...
                    'YData', [0, P(2, 2), P(3, 2)]);
    
    % Pause for animation effect
    pause(1e-3);
end

hold off;



% Functions
function dxdt = robot_dynamics(t, x, l1, l2, m1, m2, g, K, B, Q)
% Unpack state variables
q = x(1:2);
qp = x(3:4);

t1 = 4; t2 = 8; t3 = 12; t4 = 16;
   
if t<t1
    q_des = Q(1,:)';
elseif t>t1  && t<t2
    q_des = Q(2,:)';
elseif t>t2  && t<t3
    q_des = Q(3,:)';
elseif t>t3  && t<t4
    q_des = Q(4,:)';
else
    q_des = Q(5,:)';
end

% Desired velocity
qp_des = [0; 0];

% Errors
e = q - q_des;
eDot = qp - qp_des;

% Dynamics
M = mass_matrix(q(1), q(2), l1, l2, m1, m2);
G = gravity_vector(q(1), q(2), l1, l2, m1, m2, g, g);
C = coriolis_matrix(q(1), q(2), qp(1), qp(2), l1, l2, m1, m2);

% PD control with gravity compensation
Torque = -K * e - B * eDot;

% State derivatives
qpp = M \ (Torque - C * qp - G);
dxdt = [qp; qpp];
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
