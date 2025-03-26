clear; clc;
% Define desired trajectory and Middle Points
qDes = [0.1914, -0.0445, 0.3336];
xMid = [0.025, 0, 0.005; 
        0.03,  0, 0.010;
        0.035, 0, 0.015; 
        0.04,  0, 0.02; 
        0.045, 0, 0.025];
qMid = zeros(length(xMid),3);

for i = 1:length(qMid)
    qMid(i, :) = IK(xMid(i, 1), xMid(i, 2), xMid(i, 3));
end

% Parameters
time = 20;  % Total simulation time
zeta = [1 1 1];       % Third-order Prefilter Damping Ratio
wn = [1 1 1 ];        % Third-order Prefilter Natural Frequency  
kj = [50 50 50];      % Spring constants
bj = [30 30 30];      % Damping constants
wt = [0.4, 1e-5, 10]; % Weights [qDes, Time, qMid]

% Initial parameter guess
initParams = [time wn bj kj, zeta]; 

[init_T, init_Y] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj, zeta), [0 time], zeros(15, 1));
 
% Parameter bounds
lb = [5   10  1  1  10 10 10   20  20  20   0.1 0.1 0.1 ];   
ub = [10   20 20 20 40 40 40   100 100 100  3 3 3  ];  

% Objective Function
objectiveFunc = @(params) objectiveFunction(params, qDes, wt, qMid);

% Optimization setup
options = optimset('PlotFcns','optimplotfval','Display', 'off', 'TolFun', 1e-8, 'MaxIter', 400,'TolX',1e-8);
optimalParams = fmincon(objectiveFunc, initParams, [], [], [], [], lb, ub, [], options);

% Simulate with optimal parameters
[t, y] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, optimalParams(2:4), optimalParams(1), qDes, optimalParams(5:7), optimalParams(8:10),optimalParams(11:13)), [0 optimalParams(1)], zeros(15, 1));

% Plotting
disp(['Optimized Parameters: ', num2str(optimalParams)]);

% Objective function
function error = objectiveFunction(params, qDes, wt, qMid)
    x0 = zeros(15, 1);
    x0(1:3) = qDes; 

    % Simulate the system
    [t, y] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, params(2:4), params(1), qDes, params(5:7), params(8:10), params(11:13)), [0 params(1)], x0);

    % Calculate error metric
    distto1 = min(sum((y(:, 10:12) - qDes).^2, 2) + sum((params(1) - t).^2, 2)); 

    distMid = sum(arrayfun(@(i) min(sum((y(:, 10:12) - qMid(i, :)).^2, 2)), 1:size(qMid,1)));

    error = wt(1) * distto1 + wt(2) * params(1) + wt(3) * distMid;
end

% Third-order Prefiltered Dynamics
% Third-order filter modifications in the `myTwolinkwithprefilter` function:

function dxdt = myTwolinkwithprefilter(t, x, wn, time, qDes, bj, kj, zeta)
    % Define the correct size for matrix A (6x6 for the system states)
    A = [zeros(3,3), eye(3), zeros(3,3);
         -eye(3)*diag(wn)^3, -3*eye(3)*diag(zeta^2)*diag(wn)^2, -3*eye(3)*diag(zeta)*diag(wn)];
    
    % Define matrix B (6x3)
    B = [zeros(3,3); diag(wn).^3];
    
    % Extract joint positions and velocities (q = 3, qd = 3)
    q = x(7:9);     % joint positions
    qd = x(10:12);  % joint velocities
    x
    % Define control gains for spring and damping
    Kp = diag(kj);  
    Kd = diag(bj);  
    
    % Controller: PD controller with position and velocity error
    controller = Kp * (x(1:3) - q) + Kd * (x(4:6) - qd);
    
    % Get the mass, Coriolis, and gravity matrices
    [M, C, G] = compute_M_C_G(q(1), q(2), q(3), qd(1), qd(2), qd(3));
    
    % Compute joint accelerations (3x1)
    tau = M * controller + C * qd + G;
    qdd = M \ (tau - C * qd);
    
    % Return the state-space representation of the system (state derivative)
    % A * x(1:6) gives a 6x1 vector (since A is now 6x6)
    % B * qDes(:) gives a 6x1 vector
    dxdt = [A * x(1:9) + B * qDes(:); qd; qdd;qdd];  % 15x1 state vector
end



function [x, y, z] = FK(q1, q2, q3)
    l1 = 0.208; 
    l2 = 0.168;  
    x = sin(q1) .* (l1 * cos(q2) + l2 * sin(q3));
    y = l2 - l2 * cos(q3) + l1 * sin(q2);
    z = -l1 + cos(q1) .* (l1 * cos(q2) + l2 * sin(q3));
end

function Q = IK(x, y, z)
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
    Q = [q1, q2, q3] ;
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





