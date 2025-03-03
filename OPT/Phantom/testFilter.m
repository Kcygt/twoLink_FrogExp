clear; clc;

% Define desired trajectory and Middle Points
qDes = [0.1914, -0.0445, 0.3336];

% Parameters
time = 1.798620068227840;  % Time
wn = [5.31763480727411	5.14702646187068	2.35278253961713];  % Prefilter Omega     
kj = [43.0063792731132	40.4229929584824	40.2725057114223];  % Spring constants
bj = [0.108254829822687	4.99953968603787	1.70821359171333];  % Damping constants

[init_T, init_Y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj), [0 time], zeros(12, 1));

[xDes, yDes, zDes] = FK(init_Y(:,1), init_Y(:,2), init_Y(:,3));
[xAct, yAct, zAct] = FK(init_Y(:,7), init_Y(:,8), init_Y(:,9));

figure(1); hold on; grid on;
plot(xAct,zAct,'-')

plot(0.025, 0.01, '*')
plot(0.035, 0.02, '*')
plot(0.045, 0.03, '*')
% Number of time steps
num_steps = length(init_T);

% Initialize Cartesian velocity storage
cartesian_velocity = zeros(3, num_steps);

for i = 1:num_steps
    % Get joint angles at time step i
    q1 = init_Y(i,1);
    q2 = init_Y(i,2);
    q3 = init_Y(i,3);
    
    % Compute Jacobian for current joint configuration
    J = Jacobian(q1, q2, q3);
    
    % Get joint velocities at time step i
    qd = init_Y(i, 10:12)';  % Column vector of joint velocities
    
    % Compute Cartesian velocity: v = J * q_dot
    cartesian_velocity(:, i) = J * qd;
end

% Extract Cartesian velocity components
vx = cartesian_velocity(1, :);
vy = cartesian_velocity(2, :);
vz = cartesian_velocity(3, :);

% Plot Cartesian velocity over time
figure(2); hold on; grid on;
plot(init_T, vx, 'r', 'LineWidth', 1.5); hold on;
plot(init_T, vy, 'g', 'LineWidth', 1.5);
plot(init_T, vz, 'b', 'LineWidth', 1.5);
legend('Vx', 'Vy', 'Vz');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Cartesian Space Velocity');
grid on;






% Function 
function dxdt = myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj)
    zeta = 1;
    A = [zeros(3,3) eye(3);
        -eye(3)*diag(wn).^2  -eye(3)*2*zeta*diag(wn)];
    B = [zeros(3,3); diag(wn).^2];
  
    q   = x(7:9);
    qd  = x(10:12);

    Kp = -diag(kj);  
    Kd = -diag(bj);  

    [M, C, G] = compute_M_C_G(q(1), q(2), q(3), qd(1), qd(2), qd(3));

    qdd = M \ (C * qd + Kd * (qd - x(4:6)) + Kp * (q - x(1:3)));

    dxdt = [A*x(1:6) + B*qDes(:); qd; qdd];
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

function J = Jacobian(q1,q2,q3)
    % Link lenght
    l1 = 0.208; 
    l2 = 0.168;  

    % Jaccobian matrix
    J = zeros(3,3);
    J(1,1) = l1*cos(q1)*cos(q2) + l2*sin(q3)*cos(q1);
    J(1,2) = -l1*sin(q1)*sin(q2) + l2*sin(q1)*cos(q3);
    J(1,3) = l2*sin(q1)*cos(q3);

    J(2,1) = 0;
    J(2,2) = l1*cos(q2);
    J(2,3) = l2*sin(q3);

    J(3,1) = -(l1*sin(q1)*cos(q2) + l2*sin(q1)*sin(q3));
    J(3,2) = -l1*sin(q2)*cos(q1);
    J(3,3) = l2*cos(q1)*cos(q3);

end


function dvar=quickdiff(tt,var,method)
% diff_of_var=QUICKDIFF(timevec,var,method)
% A quick way to calculate the differential. 
% Estimated as $dx/dt \approx \frac{x_{n+1}-x_n}{t_{n+1}-t_n)$ and
% this is interpolated back on to the time vector. Thus the time vector
% does not need to have equal increments.
%
% Best to convert the time vector to numbers before the call
%
% Most of the code is handling time structures, the actual work is
% done with midpts=... and dvar=... This could be simplified!
%
% 'method' is not implimented for now, but could just be the method
% passed to interp1
%
% can var be a matrix?

    if isduration(tt) || isdatetime(tt)
        warning('Best to convert the time vector to numbers');
    end

    deltaTimes=diff(tt);
    if isduration(deltaTimes)
        if ~strcmp(deltaTimes.Format,'s') % will assume seconds
            warning(sprintf('Time vector format is %s',deltaTimes.Format))
        end
        if isdatetime(tt)
            tt=timeofday(tt);
        end
        tt=seconds(tt); % finally a number?        
        midpts=(tt(1:end-1)+tt(2:end))/2; % find midpoints in time vector
        dvar=interp1(midpts,diff(var)./seconds(deltaTimes),tt); % interpolate diffs to the original time vector
    else % assume time is just numbers
        midpts=(tt(1:end-1)+tt(2:end))/2; % find midpoints in time vector

        dvar=interp1(midpts,diff(var)./diff(tt),tt); % interpolate diffs to the original time vector
    end
end

