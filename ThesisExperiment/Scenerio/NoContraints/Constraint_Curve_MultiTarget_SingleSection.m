clear; clc;
close all;

% Define desired trajectory and Middle Points
qDes = [0.1914, -0.0445, 0.3336];
[xDes, yDes, zDes] = FK(qDes(1), qDes(2), qDes(3));
xDes = [xDes, yDes, zDes];

xMid = zeros(2,3);
qMid = zeros(2,3);

xMid(1,:) = [0.03, 0, 0.01];
xMid(2,:) = [0.04, 0, 0.02];
qMid(1,:) = IK(xMid(1,1), xMid(1,2), xMid(1,3));
qMid(2,:) = IK(xMid(2,1), xMid(2,2), xMid(2,3));

% Parameters
tspan = 10;
zeta = [1 1 1];
wn = [1 1 1];
Kp = [50 50 50];
Kd = [20 20 20];

% Weights
wt = [350, 0.01, 0.002]; % [Target, End, Time]

initPrms = [tspan, zeta,wn,Kp,Kd];

% Initial Condition
[ti, yi] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, tspan, qDes,zeta, wn,Kp,Kd), [0 tspan], zeros(12, 1));

% Lower and Upper Limits
lb = [2 ... % time
      0.4 0.4 0.4 ... % zeta
      0.5 0.5 0.5 ... % Wn
      30 30 30 ...    % Kp
      10 10 10];      % Kd 
ub = [5 ... % time
      1.1 1.1 1.1 ... % zeta
      50 50 50 ...    % Wn
      200 200 200 ... % Kp
      70 70 70];      % Kd

% Objective Function
objectiveFunc = @(params) objectiveFunction(params, qDes, wt, xMid, xDes);

% Run optimization
options = optimset('PlotFcns', 'optimplotfval', 'Display', 'off', 'TolCon', 1e-6); % Added constraint tolerance

[Opt, fval] = fmincon(objectiveFunc, initPrms, [], [], [], [], lb, ub, ...
                      @(prms) trajConstraint(prms, qDes, xMid), options);

% Simulate with optimal parameters
[tt, yy] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, Opt(1), qDes, Opt(2:4),Opt(5:7),Opt(8:10),Opt(11:13)), ...
                  [0 Opt(1)], zeros(12, 1));

%%% Plotting
[xi, yi_plot, zi] = FK(yi(:,7), yi(:,8), yi(:,9)); % Initial Trajectory
[x_opt, y_opt, z_opt] = FK(yy(:,7), yy(:,8), yy(:,9)); % Optimized Trajectory

figure; hold on; grid on;
plot(xi, zi,'--')
plot(x_opt,z_opt,'.-')
plot(xMid(1,1),xMid(1,3),'*')
plot(xMid(2,1),xMid(2,3),'*')

plot(xDes(1),xDes(3),'o')
legend('Initial Trajectory','Optimized Trajectory','Midpoint','Endpoint')

disp('Optimal Parameter:')
disp(['Time: ', num2str(Opt(1))])
disp(['Zeta: ', num2str(Opt(2:4))])
disp(['Wn: ', num2str(Opt(5:7))])
disp(['Kp: ', num2str(Opt(8:10))])
disp(['Kd: ', num2str(Opt(11:13))])

% Objective Function
function error = objectiveFunction(prms, qDes, wt, xMid, xDes)
    x0 = zeros(12, 1);
    x0(1:3) = qDes;

    % Simulate the system
    [~, y] = ode23s(@(t,x) myTwolinkwithprefilter(t,x,prms(1),qDes,prms(2:4),prms(5:7),prms(8:10),prms(11:13)), ...
                    [0 prms(1)], x0);

    [xOut,~,zOut] = FK(y(:,7),y(:,8),y(:,9));
    
    % Calculate minimum distance to middle point
    dx1 = abs(xOut - xMid(1,1)).^2;
    dz1 = abs(zOut - xMid(1,3)).^2;
    distMid1 = sqrt(dx1+dz1);
    
    dx2 = abs(xOut - xMid(2,1)).^2;
    dz2 = abs(zOut - xMid(2,3)).^2;
    distMid2 = sqrt(dx2+dz2);

    % End point error
    dxEnd = abs(xOut(end) - xDes(1)).^2;
    dzEnd = abs(zOut(end) - xDes(3)).^2;
    distEndErr = sqrt(dxEnd + dzEnd);
    
    % Time penalty
    timePenalty = prms(1);

    % Composite error (normalized)
    error = wt(1)*min(distMid1+distMid2) + wt(2)*distEndErr + wt(3)*timePenalty;
end

% Constraint Function for Midpoint Proximity
function [c, ceq] = trajConstraint(prms,qDes,xMid)
    ceq = []; % No equality constraints

    % Simulate trajectory
    [~, yy] = ode23s(@(t,x) myTwolinkwithprefilter(t,x,prms(1),qDes,prms(2:4),prms(5:7),prms(8:10),prms(11:13)), ...
                    [0 prms(1)], zeros(12, 1));
    [x,y,z] = FK(yy(:,7),yy(:,8),yy(:,9));     % Optimized Trajectory
    
    % Calculate distances to midpoint in 3D space
    dx1 = (abs(x - xMid(1,1)).^2 + abs(z - xMid(1,3)).^2);
    dx2 = (abs(x - xMid(2,1)).^2 + abs(z - xMid(2,3)).^2);

    distance1 = sqrt(dx1);
    distance2 = sqrt(dx2);
    % End point error
    dxEnd = abs(x(end) - 0.05).^2;
    dzEnd = abs(z(end) - 0.05).^2;
    distEndErr = sqrt(dxEnd + dzEnd);
    
    % Nonlinear inequality constraint: min distance <= 10cm (0.1m)
    c = [min(distance1) - 0.0005;
         min(distance2) - 0.0005;
         distEndErr - 0.0005]; 
end

% Dynamics Function with Prefilter
function dxdt= myTwolinkwithprefilter(t,x,t_st,qDes,zeta,wn,Kp,Kd)
    % Filter 
    A1=[zeros(3), eye(3); -diag(wn).^2,-2*diag(zeta)*diag(wn)];
    B1=[zeros(3); diag(wn).^2];
    % Actual Joint position and velocity
    q=x(7:9);
    qd=x(10:12);
    % Controller
    controller=diag(Kp) * (x(1:3)-q) + ...
               diag(Kd) * (x(4:6)-qd);
    % Robot Dynamics
    [M,C,G]=compute_M_C_G(q(1),q(2),q(3),qd(1),qd(2),qd(3));
    % Torque
    tau=M*(controller)+C*qd;
    % Acceleration
    qdd=M\(tau-C*qd);
    % Output
    dxdt=[A1*x(1:6)+B1*qDes(:); qd; qdd];
end

% Forward Kinematics (FK)
function [x,y,z]=FK(q1,q2,q3)
    l1=0.208; 
    l2=0.168;  
    x=sin(q1).*(l1*cos(q2)+l2*sin(q3));
    y=l2-l2*cos(q3)+l1*sin(q2);
    z=-l1+cos(q1).*(l1*cos(q2)+l2*sin(q3));
end

% Inverse Kinematics (IK)
function Q=IK(x,y,z)
    l1=0.208; 
    l2=0.168;  
    q1=atan2(x,z+l1);

    R=sqrt(x^2+(z+l1)^2);
    r=sqrt(x^2+(y-l2)^2+(z+l1)^2);
    
    Beta=atan2(y-l2,R);
    Gamma=acos((l1^2+r^2-l2^2)/(2*l1*r));
    
    q2=Gamma+Beta;

    Alpha=acos((l1^2+l2^2-r^2)/(2*l1*l2));
    
    q3=q2+Alpha-pi/2;
    
    Q=[q1,q2,q3];
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


