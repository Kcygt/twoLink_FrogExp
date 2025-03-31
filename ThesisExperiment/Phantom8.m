clear; clc;
% close all;
% Define desired trajectory and Middle Points
qDes = [0.1914, -0.0445, 0.3336];
xMid = zeros(4,3);
xMid(1,:) = [0.01, 0, 0.01 ];
xMid(2,:) = [0.02, 0, 0.03 ];
xMid(3,:) = [0.0365, 0, 0.034];
xMid(4,:) = [0.04, 0, 0.045 ];

qMid = zeros(4,3);
qMid(1,:) = IK(xMid(1,1), xMid(1,2), xMid(1,3));
qMid(2,:) = IK(xMid(2,1), xMid(2,2), xMid(2,3));
qMid(3,:) = IK(xMid(3,1), xMid(3,2), xMid(3,3));
qMid(4,:) = IK(xMid(4,1), xMid(4,2), xMid(4,3));


% Parameters
num_stages = length(qMid) + 1; % Number of different parameter sets

time_stages = [0.85, 1.35, 3.85,4.95,6 ];
tspan = [0, time_stages(end)];

zeta = [1 1 1; 0.2 1 0.5;    .6 1 0.8;     0.01 1 .01;  1 1 1]; 
wn   = [1 1 1; 1   1  2;     10 1 0.2;     1 1  0.5;   1 1 1];
kj = [60 50 40; 50 50 50; 50 50 50;50 50 50;50 50 50];
bj = [30 30 30; 30 30 30; 30 30 30;30 30 30;30 30 30];

wt = [0.5, 1e-5, 200];  % Weights [qDes, Time, qMid]

[init_T, init_Y] = ode45(@(t, x) myTwolinkwithprefilter(t, x, wn, zeta, time_stages, qDes, bj, kj), tspan,  zeros(12, 1));

%%% Plotting
[x,y,z] = FK(init_Y(:,7),init_Y(:,8),init_Y(:,9));
figure(1); hold on; grid on;
plot(x,z)
plot(xMid(1,1),xMid(1,3),'*')
plot(xMid(2,1),xMid(2,3),'*')
plot(xMid(3,1),xMid(3,3),'*')
plot(xMid(4,1),xMid(4,3),'*')

plot(0.05,0.05,'o')
%%%%
% 
% figure(2); hold on; grid on;
% plot(init_T,init_Y(:,1))
% plot(init_T,init_Y(:,3))

% myTwolinkwithprefilter function
function dxdt= myTwolinkwithprefilter(t, x, wn, zeta, t_st, qDes, bj, kj)
    for i = 1:length(wn)
        W = diag(wn(i,:));
        Z = diag(zeta(i,:));
        Kp{i} = diag(kj(i,:));
        Kd{i} = diag(bj(i,:));
        q = x(7:9);
        qd = x(10:12);
    
        A{i} = [zeros(3), eye(3); -W.^2, -2*Z*W];
        B{i} = [zeros(3); W.^2];
        
        controller{i} = Kp{i} * (x(1:3) - q) + Kd{i} * (x(4:6) - qd);
        [M, C, ~] = compute_M_C_G(q(1), q(2), q(3), qd(1), qd(2), qd(3));
        
        tau{i} = M * controller{i} + C * qd;
        qdd{i} = M \ (tau{i} - C * qd);
        output{i} = [A{i} * x(1:6) + B{i} * qDes(:); qd; qdd{i}];
    end


    if t < t_st(1)
        dxdt = output{1};
    elseif t_st(1) <= t && t < t_st(2)
        dxdt = output{2};
    elseif t_st(2) <= t && t < t_st(3)
        dxdt = output{3};
    elseif t_st(3) <= t && t < t_st(4)
        dxdt = output{4};
    else
        dxdt = output{5};
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





