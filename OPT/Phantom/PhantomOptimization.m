% Optimization of two-link robot arm tracking
clear; clc;

% Define desired trajectory and Middle Points
% qDes = [ 0.4481, -0.0075, 0.1361];
qDes = [ 0.7854, pi/4, pi/8 ];
qMid = zeros(4,3);
qMid(1,:) = IK(0.03, 0,-0.0127);
qMid(2,:) = IK(0.06, 0,-0.0268);
qMid(3,:) = IK(0.09, 0,-0.0409);
qMid(4,:) = IK(0.12, 0,-0.0549);

%  Parameters  0.0528404      10.2888      11.2778      10.7457      0.19054      5.13513      9.99989      8.29814      50.5053      99.9982
time = 10;      % time
wn = [ 5 3 5];              % Prefilter Omega     
kj = [10 10 10];        % Spring  [q1 q2]
bj = [0.3 0.3  0.9];        % Damping [q1 q2]
wt = [10, .1, 18000];         % weights [qDes, Time, qMid]

% Optimization setup
initParams = [time  wn bj kj]; % Initial guess for [time, wn, bj, kj]

[init_T, init_Y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj),   [0 time], zeros(12, 1));
[xInit,yInit,zInit] = FK(init_Y(:,7),init_Y(:,8),init_Y(:,9));
plot3(xInit,yInit,zInit,'.')
% figure(1);
% plot(init_Y(:,1),'*-')
% 
% figure(2);
% plot(init_Y(:,3),'*-')
% 
% figure(3); hold on;
% plot(xInit,zInit,'-.')
% plot(0.12, -0.0549,'o')
% plot(0.09, -0.0409,'o')
% plot(0.06, -0.0268,'o')
% plot(0.03, -0.0127,'o')

% Lower and upper boundaries 
lb = [1     0.01   0.01  0.01    .001  .001  .001  0.01  0.01 0.01];   % Lower bounds
ub = [2     30    30 30     5 5 5      10 10 10 ];     % Upper bounds

% Objective Function
objectiveFunc = @(params) objectiveFunction(params, qDes, wt, qMid);

% Run optimization
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'MaxIter', 400);
optimalParams = fmincon(objectiveFunc, initParams, [], [], [], [], lb, ub, [], options);

% Simulate with optimal parameters and plot results
[t, y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, optimalParams(2:4), optimalParams(1), qDes, optimalParams(5:7), optimalParams(8:10)), [0 optimalParams(1)], zeros(12, 1));

% Output
[xAct,yAct,zAct] = FK(y(:,7),y(:,8),y(:,9));
[xDes,yDes,zDes] = FK(qDes(:, 1), qDes(:, 2), qDes(:,3));
% [xInit,yInit,zInit] = FK(init_Y(:,7),init_Y(:,8),init_Y(:,9));

% Plotting
% Desired, Actual and Optimised Data
figure(1); hold on; grid on;
plot(xInit, zInit, '-');
legend('Initial')
title('Initial Trajectory Tracking');

figure(2); hold on; grid on;
plot(xAct, zAct, '-');
plot(xDes, zDes, 'o-');
plot(0.03, -0.0127,'*')
plot(0.06, -0.0268,'*')
plot(0.09, -0.0409,'*')
plot(0.12, -0.0549,'*')

xlabel('X axis'); ylabel('Y axis');
legend('Optimised', 'Desired');
title('Optimized Trajectory Tracking');
disp(['Optimized Parameters :', num2str(optimalParams)])

% % Mid points in joint space
% figure(2);plot(y(:,5),y(:,6),qMid(1,:),qMid(2,:),'o');
% xlabel('Joint 1 position')
% ylabel('Joint 2 position')
% 
% title('joint space of a (near) optimal staight line in cartesian space')
% 
% % joint space plot
% figure(3); grid on; hold on;
% plot(t,y(:,5:6));
% xlabel('Time (s)')
% ylabel('Position (rad)')
% legend('Q1','Q2')
% title('Joint position (rad)')
% 
% %  cartesian space plot
% figure(4); hold on; grid on;
% plot(xAct(:,1),xAct(:,2))
% xlabel('X axis')
% ylabel('Y axis')
% legend('X','Y')
% title('Cartesian Position (m)')
% 
% % x/y vs time
% figure(5); grid on; hold on;
% plot(t,xAct(:,1:2))
% xlabel('Time (s)')
% ylabel('Position')
% legend('X','Y')
% title('Cartesian Position vs Time')

% publish('simOpt.m','pdf');
% disp(sprintf('KY %s \t %s \t %s',mfilename,pwd,datetime("now")));

% Objective function
function error = objectiveFunction(params, qDes,wt,qMid)
    
    % Initial conditions
    x0 = zeros(12, 1);
    x0(1:3) = [qDes(1, 1); qDes(1, 2);qDes(1, 3)];
    
    % Simulate the system
    [t, y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, params(2:4), params(1), qDes, params(5:7), params(8:10)),   [0 params(1)], x0);

    % Calculate the error metric 
    distto1 = min(sum((y(:, 7:9) - qDes(1,:)).^2,2) + sum((params(1) - t).^2,2)); 

    distMid1 = min(sum((y(:, 7:9) - qMid(1,:)).^2,2));  
    distMid2 = min(sum((y(:, 7:9) - qMid(2,:)).^2,2));       
    distMid3 = min(sum((y(:, 7:9) - qMid(3,:)).^2,2));       
    distMid4 = min(sum((y(:, 7:9) - qMid(4,:)).^2,2));      

    % time1 = min(sum((params(1) - t).^2,2));
    % time2 = min(sum((params(2) - t).^2,2));
    time1 = params(1);

    error   = wt(1) * distto1  + ...  % Desired
              wt(2) * time1    + ...  % time
              wt(3) * distMid1 + wt(3) * distMid2 + ...  % Mid-point
              wt(3) * distMid3 + wt(3) * distMid4;  % Mid-point


    % distto5 = 5000 * sum((y(:, 5:6) - qMid3'),2) + w2 * (sum(  (   (time(1) + (time(2) - time(1))/2 ) - t).^2   ,2));

end

% myTwolinkwithprefilter function
function dxdt = myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj)
    zeta = 1;
    A = [zeros([3 3]) eye(3); -eye(3)*wn(1)^2 -eye(3)*2*zeta*wn(1)];
    B = [0 0 0; 0 0 0;0 0 0; wn(1)^2 0 0;0 wn(2)^2 0; 0 0 wn(3)^2 ];
    

    % Actual position and velocity
    q   = x(7:9);
    qd  = x(10:12);
    
    q1p = x(10); q2p = x(11); q3p = x(12);
    q1  = x(7); q2 = x(8); q3 = x(9);
    
    % Robot constants
    [M,C,G] = compute_M_C_G(q1,q2,q3,q1p,q2p,q3p);
    Numerator = C*qd + [-bj(1) 0 0; 0 -bj(2) 0; 0 0 -bj(3)  ]*(qd-x(4:6)) + [-kj(1) 0 0; 0 -kj(2) 0; 0 0 -kj(3)]*(q - x(1:3));
    qdd = M\Numerator;
   
    dotx = A*x(1:6) + B*qDes(1, :)';
   
    dxdt = [dotx; qd; qdd];
end

function [x,y,z] = FK(q1,q2,q3)
    l1 = 0.208; 
    l2 = 0.168;  
    
    x = sin(q1) .* (l1*cos(q2) + l2*sin(q3));
    y = l2-l2*cos(q3)+l1*sin(q2);
    z = -l1+cos(q1).*(l1*cos(q2)+l2*sin(q3));
end

function [q1,q2,q3] = IK(x,y,z)
    l1 = 0.208; 
    l2 = 0.168;  
    q1 = atan2(x,z+l1);
    
    R = sqrt(x^2 + (z+l1)^2);
    r = sqrt(x^2 + (y-l2)^2 + (z+l1)^2);
    Beta  = atan2(y-l2,R);
    Gamma =  acos((l1^2+r^2 - l2^2)/(2*l1*r));
    q2 = Gamma + Beta;

    Alpha = acos((l1^2 + l2^2 - r^2)/ (2*l1*l2));
    q3 = q2 + Alpha - pi/2;
end


function [M, C, G] = compute_M_C_G(theta1, theta2,theta3, dtheta1, dtheta2,dtheta3)
    % link lenghts
    l1 = 0.208;   l2 = 0.168;   l3 = 0.0325;
    l5 = -0.0368; l6 = 0.0527;
    
    % Gravity
    g = -9.80665; % m/s^2
       
    % segment A
    m_a = 0.0202;
    Ia_xx = 0.4864*1e-4;  Ia_yy = 0.001843*1e-4;  Ia_zz = 0.4864*1e-4;
    Ia = [Ia_xx 0 0;  0  Ia_yy 0;  0 0 Ia_zz];
    
    % segment C
    m_c = 0.0249;
    Ic_xx = 0.959*1e-4;  Ic_yy = 0.959*1e-4;  Ic_zz = 0.0051*1e-4;
    Ic = [Ic_xx 0 0;  0  Ic_yy 0;  0 0 Ic_zz];
    
    % segment BE
    m_be = 0.2359;
    Ibe_xx = 11.09*1e-4;  Ibe_yy = 10.06*1e-4;  Ibe_zz = 0.591*1e-4;
    Ibe = [Ibe_xx 0 0;  0  Ibe_yy 0;  0 0 Ibe_zz];
    
    % segment DF
    m_df = 0.1906;
    Idf_xx = 7.11*1e-4;  Idf_yy = 0.629*1e-4;  Idf_zz = 6.246*1e-4;
    Idf = [Idf_xx 0 0;  0  Idf_yy 0;  0 0 Idf_zz];
    
    % BASE
    Ibaseyy = 11.87e-4;
    
    % MASS 
    M11 = ( 1/8*( 4*Ia_yy  + 4*Ia_zz  + 8*Ibaseyy + 4*Ibe_yy + 4*Ibe_zz + 4*Ic_yy + 4*Ic_zz + 4*Idf_zz + 4*l1^2*m_a + l2^2*m_a + l1^2*m_c + 4*l3^2*m_c  ) + ...
            1/8*( 4*Ibe_yy - 4*Ibe_zz + 4*Ic_zz   + l1^2*(4*m_a + m_c)) * cos(2*theta2) + ...
            1/8*( 4*Ia_yy  - 4*Ia_zz  + 4*Idf_yy  - 4*Idf_zz - l2^2*m_a - 4*l3^2*m_c) * cos(2*theta3) + l1*(l2*m_a + l3*m_c)*cos(theta2)*sin(theta3)  );
    
    M22 = 1/4*(4*(Ibe_xx + Ic_xx + l1^2*m_a) + l1^2*m_c);
    M23 = -1/2*l1*(l2*m_a + l3*m_c) * sin(theta2-theta3);
    M32 = M23;
    M33 = 1/4 * (4*Ia_xx + 4*Idf_xx + l2^2*m_a + 4*l3^2*m_c);
    
    M = [M11 0 0; 0 M22 M23; 0 M32 M33];
    
    % CORIOLIS
    C11 = 1/8*( -2*sin(theta2) * ((4*Ibe_yy - 4*Ibe_zz * 4*Ic_yy - 4*Ic_zz + 4*l1^2*m_a +l1^2*m_c )*cos(theta2) + 2*l1*(l2*m_a + l3*m_c)*sin(theta3) ) )*dtheta2 + ...
           2*cos(theta3)*(2*l1*(l2*m_a + l3*m_c)*cos(theta2) + (-4*Ia_yy + 4*Ia_zz -4*Idf_yy + 4*Idf_zz + l2^2*m_a + 4*l3^2*m_c)*sin(theta3)) * dtheta3 ;
    
    C12 = -1/8 * ((4*Ibe_yy - 4*Ibe_zz + 4*Ic_yy - 4*Ic_zz + l1^2*(4*m_a + m_c))*sin(2*theta2) + 4*l1*(l2*m_a+l3*m_c)*sin(theta2)*sin(theta3))*dtheta1; 
    C13 = -1/8*(-4*l1*(l2*m_a + l3*m_c)*cos(theta2)*cos(theta3) - (-4*Ia_yy + 4*Ia_zz - 4*Idf_yy + 4*Idf_zz + l2^2*m_a + 4*l3*m_c)*sin(2*theta3))*dtheta1;
    C21 = -C12;
    C23 = 1/2*l1*(l2*m_a + l3*m_c)*cos(theta2*theta3)*dtheta3;
    C31 = -C13;
    C32  = -1/2*l1*(l2*m_a + l3*m_c)*cos(theta2 - theta3)*dtheta2;
    
    C = [C11 C12 C13; C21 0 C23; C31 C32 0];
    
    % GRAVITY
    N2 = 1/2*g*(2*l1*m_a + 2*l5*m_be + l1*m_c)*cos(theta2);
    N3 = 1/2*g*(l2*m_a + 2*l3*m_c - 2*l6*m_df)*sin(theta3);
    
    G = [0 N2 N3]';
end

