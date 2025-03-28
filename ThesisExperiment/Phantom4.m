clear; clc;
% close all;
% Define desired trajectory and Middle Points
qDes = [0.1914, -0.0445, 0.3336];
xMid = zeros(2,3);
xMid(1,:) = [0.01, 0, 0.01 ];
xMid(2,:) = [0.03, 0, 0.06 ];

qMid = zeros(2,3);
qMid(1,:) = IK(xMid(1,1), xMid(1,2), xMid(1,3));
qMid(2,:) = IK(xMid(2,1), xMid(2,2), xMid(2,3));



% Parameters
time = 30;  % Time
t1 = 0.7;
t2 = 6;
t3 = 9;
zeta1 = [1 1 1];       % Prefilter Zeta
zeta2 = [.3 1 .5];       % Prefilter Zeta
zeta3 = [.9 1 .9];       % Prefilter Zeta

wn1 = [1 1 1 ];          % Prefilter Omega     
wn2 = [11 1 1 ];          % Prefilter Omega     
wn3 = [12 1 1 ];          % Prefilter Omega     

kj1 = [60 50 40];       % Spring constants
bj1 = [30 30 30];       % Damping constants
kj2 = [50 50 50];       % Spring constants
bj2 = [30 30 30];       % Damping constants
kj3 = [50 50 50];       % Spring constants
bj3 = [30 30 30];       % Damping constants

wt = [0.5, 1e-5, 200];  % Weights [qDes, Time, qMid]

% Optimization setup
initParams = [t1 t2 t3 wn1 wn2 wn3 bj1 bj2 bj3 kj1 kj2 kj3 zeta1 zeta2 zeta3]; % Initial guess

[init_T, init_Y] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, wn1,wn2,wn3 ,t1,t2,t3, qDes, bj1,bj2,bj3, kj1,kj2,kj3, zeta1,zeta2,zeta3), [0 t3], zeros(12, 1));

%%% Plotting
% [x,y,z] = FK(init_Y(:,7),init_Y(:,8),init_Y(:,9));
% figure(1); hold on; grid on;
% plot(x,z)
% plot(xMid(1,1),xMid(1,3),'*')
% plot(xMid(2,1),xMid(2,3),'*')

%%%%
% Upper and Lower Limits
lb = [0  0  0     10  1  1   11  1  1    12 1 1          10 15 16  10 10 10 10 10 10        20  20  20  20  20  20  20  20  20       0 0 0 0 0 0 0 0 0 ];   
ub = [10 10 10    20 20 20   20 20 20    20 20 20        40 40 40  40 40 40 40 40 40        100 100 100 100 100 100 100 100 100      1 1 1 1 1 1 1 1 1];  

% Objective Function
objectiveFunc = @(params) objectiveFunction(params, qDes, wt, qMid);

% Run optimization
options = optimset('PlotFcns','optimplotfval','Display', 'off', 'TolFun', 1e-8, 'MaxIter', 400,'TolX',1e-8);
optimalParams = fmincon(objectiveFunc, initParams, [], [], [], [], lb, ub, [], options);

% Simulate with optimal parameters
[t, y] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, optimalParams(4:6),optimalParams(7:9),optimalParams(10:12), optimalParams(1),optimalParams(2),optimalParams(3), ...
    qDes, optimalParams(13:15),optimalParams(16:18),optimalParams(19:21),optimalParams(22:24),optimalParams(25:27),optimalParams(28:30)...
    ,optimalParams(31:33),optimalParams(34:36),optimalParams(37:39)), [0 optimalParams(3)], zeros(12, 1));

Plotting
disp(['Optimized Parameters: ', num2str(optimalParams)]);

% Objective function
function error = objectiveFunction(params, qDes, wt, qMid)
    x0 = zeros(12, 1);
    x0(1:3) = qDes; 

    % Simulate the system
    [t, y] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, params(4:6),params(7:9),params(10:12), params(1),params(2),params(3), ...
    qDes, params(13:15),params(16:18),params(19:21),params(22:24),params(25:27),params(28:30)...
    ,params(31:33),params(34:36),params(37:39)), [0 params(2)], x0);

    % Calculate error metric
    distto1 = min(sum((y(:, 7:9) - qDes).^2, 2) + sum((params(1) - t).^2, 2)); 

    distMid = sum(arrayfun(@(i) min(sum((y(:, 7:9) - qMid(i, :)).^2, 2)), 1:size(qMid,1)));

    error = wt(1) * distto1 + wt(2) * params(3) + wt(3) * distMid;

end

% myTwolinkwithprefilter function
function dxdt= myTwolinkwithprefilter(t, x, wn1,wn2,wn3, t1,t2, t3, qDes,  bj1,bj2,bj3,   kj1,kj2,kj3,   zeta1,zeta2,zeta3)
    % zeta = 1;

    A1 = [zeros(3,3) eye(3);
        -eye(3)*diag(wn1).^2  -eye(3)*2*diag(zeta1)*diag(wn1)];
    B1 = [zeros(3,3); diag(wn1).^2];

    A2 = [zeros(3,3) eye(3);
        -eye(3)*diag(wn2).^2  -eye(3)*2*diag(zeta2)*diag(wn2)];
    B2 = [zeros(3,3); diag(wn2).^2];

    A3 = [zeros(3,3) eye(3);
        -eye(3)*diag(wn3).^2  -eye(3)*2*diag(zeta3)*diag(wn3)];
    B3 = [zeros(3,3); diag(wn3).^2];
    

    q   = x(7:9);
    qd  = x(10:12);
    Kp1 = diag(kj1);  
    Kd1 = diag(bj1);  
    
    Kp2 = diag(kj2);  
    Kd2 = diag(bj2);  

    Kp3 = diag(kj3);  
    Kd3 = diag(bj3);  
    
    controller1 = Kp1 * (x(1:3) - q) + Kd1 * (x(4:6) - qd);
    controller2 = Kp2 * (x(1:3) - q) + Kd2 * (x(4:6) - qd);
    controller3 = Kp3 * (x(1:3) - q) + Kd3 * (x(4:6) - qd);

    [M, C, G] = compute_M_C_G(q(1), q(2), q(3), qd(1), qd(2), qd(3));

    tau1 = M * (controller1) + C * qd ;
    tau2 = M * (controller2) + C * qd ;
    tau3 = M * (controller3) + C * qd ;

    qdd1 = M \ (tau1 - C * qd );
    qdd2 = M \ (tau2 - C * qd );
    qdd3 = M \ (tau3 - C * qd );

    % qdd = M \ ( torque - C * qd);
    output1 = [A1*x(1:6) + B1*qDes(:); qd; qdd1];
    output2 = [A2*x(1:6) + B2*qDes(:); qd; qdd2];
    output3 = [A3*x(1:6) + B3*qDes(:); qd; qdd3];
    
    if t < t1
        dxdt = output1;
    elseif t1 < t && t < t2
        dxdt = output2;
    else
        dxdt = output3;
    end


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





