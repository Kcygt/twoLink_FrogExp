clear; clc;
close all;

% Define desired trajectory and Middle Points
qDes = [0.1914, -0.0445, 0.3336];
[xDes, yDes, zDes] = FK(qDes(1),qDes(2),qDes(3));
xDes = [xDes, yDes, zDes];
xMid = zeros(2,3);
xMid(1,:) = [0.02, 0, 0.01 ];
xMid(2,:) = [0.04, 0, 0.03 ];
qMid = zeros(2,3);
qMid(1,:) = IK(xMid(1,1), xMid(1,2), xMid(1,3));
qMid(2,:) = IK(xMid(2,1), xMid(2,2), xMid(2,3));

% Parameters
tspan = [10 20 30];
zeta1 = [1 1 1];
zeta2 = [1 1 1];
zeta3 = [1 1 1];
wn1 = [1 1 1]; 
wn2 = [1 1 1];
wn3 = [1 1 1];

% weights
wt = [250, 1e+5, .5];  % [Target, End, Time]
% wt = [150, 1e+5, .001];  % [Target, End, Time]

initPrms = [tspan,zeta1,zeta2,zeta3,wn1,wn2,wn3];

% Initial Condition
[ti, yi] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, tspan, qDes, zeta1, zeta2,zeta3, wn1,wn2,wn3),[0 tspan(3)], zeros(12, 1));


% Lower and Upper Limits
lb = [1 1 1  ...         % time 
      0.1 0.1 0.1  ...    % zeta1
      0.1 0.1 0.1  ...    % zeta2
      0.1 0.1 0.1  ...    % zeta3
      0.1 0.1 0.1  ...    % Wn1
      0.1 0.1 0.1 ...     % wn2
      0.1 0.1 0.1  ];     % Wn3
ub = [10 10 10 ...            % time
      1 1 1 ...       % zeta1
      1 1 1 ...       % zeta2
      1 1 1 ...       % zeta2
      5 5 5 ...    % wn1
      5 5 5 ...    % wn2
      5 5 5];      % wn3


% Objective Function
objectiveFunc = @(params) objectiveFunction(params, qDes, wt, qMid,xMid,xDes);


% Run optimization
options = optimset('PlotFcns', 'optimplotfval', 'Display', 'off');
% [Opt,fval] = fmincon(objectiveFunc, initPrms, [], [], [], [], lb, ub, ...
%                     @(prms) trajConstraint(prms, qDes, xMid), options);

[Opt,fval] = fmincon(objectiveFunc, initPrms, [], [], [], [], lb, ub, ...
                    @(prms) trajConstraint(prms, qDes, xMid, xDes), options);


% Simulate with optimal parameters
[tt, yy] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, Opt(1:3), qDes, Opt(4:6),   Opt(7:9),Opt(10:12),Opt(13:15),Opt(16:18),Opt(19:21)),   ...
                                                                    [0 Opt(3)], zeros(12, 1));

%%% Plotting
[xi,yi,zi] = FK(yi(:,7),yi(:,8),yi(:,9));  % Initial Trajectory
[x,y,z] = FK(yy(:,7),yy(:,8),yy(:,9));     % Optimized Trajectory

figure; hold on; grid on;
plot(xi,zi,'--')
plot(x,z,'.-')
plot(xMid(1,1), xMid(1,3),'*') 
plot(xMid(2,1), xMid(2,3),'*') 
plot(0.05,0.05,'o')
legend('Initial Trajectory','Optimized Trajectory')

disp('Optimal Parameter:')
disp(['Time: ', num2str(Opt(1:3))])
disp(['Zeta1: ', num2str(Opt(4:6))])
disp(['Zeta2: ', num2str(Opt(7:9))])
disp(['Zeta3: ', num2str(Opt(10:12))])
disp(['Wn1: ', num2str(Opt(13:15))])
disp(['Wn2: ', num2str(Opt(16:18))])
disp(['Wn3: ', num2str(Opt(19:21))])


function [c, ceq] = trajConstraint(prms, qDes, xMid, xDes)
    % Simulate trajectory
    [t, y] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, prms(1:3), qDes, prms(4:6), prms(7:9),prms(10:12),prms(13:15),prms(16:18),prms(19:21)), ...
                    [0 prms(3)], zeros(12, 1));
    [xOut, yOut, zOut] = FK(y(:,7), y(:,8), y(:,9));
    
    % Middle point constraint
    distances1 = sqrt( (xOut - xMid(1,1)).^2 + (yOut - xMid(1,2)).^2 + (zOut - xMid(1,3)).^2);
    distances2 = sqrt( (xOut - xMid(2,1)).^2 + (yOut - xMid(2,2)).^2 + (zOut - xMid(2,3)).^2);
    minDist1 = min(distances1);
    minDist2 = min(distances2);
    
    % Endpoint constraint
    finalPos = [xOut(end) yOut(end) zOut(end)];
    endError = norm(finalPos - xDes);
    
    % Combined inequality constraints
    c = [minDist1 - 0.005;
         minDist2 - 0.005;
         endError - 0.005;
         max(y(:,10)) - 0.01;
         max(y(:,11)) - 0.01;
         max(y(:,12)) - 0.01;
         prms(1) - prms(2);
         prms(2) - prms(3)];   % Final position error < 10cm
    ceq = [];
end


function error = objectiveFunction(prms, qDes, wt, qMid, xMid, xDes)
    x0 = zeros(12, 1);
    x0(1:3) = qDes; 
    
    % Simulate the system
    [t, y] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, prms(1:3), qDes, prms(4:6), prms(7:9),prms(10:12),prms(13:15),prms(16:18),prms(19:21)), ...
                    [0 prms(3)], x0);
    
    [xOut, yOut, zOut] = FK(y(:,7), y(:,8), y(:,9));
    xOut = [xOut, yOut, zOut];
    
    % Calculate minimum distance to middle point
    distances1 = sqrt(sum((xOut - xMid(1,:)).^2, 2));
    distances2 = sqrt(sum((xOut - xMid(2,:)).^2, 2));

    [minDist1, idx] = min(distances1) ;

    [minDist2, idx] = min(distances2) ;
    
    
    % End point error
    endError = norm(xOut(end,:) - xDes);
    
    % Time penalty
    timePenalty = prms(3);
    
    % Composite error
    error = wt(1) * (minDist1 + minDist2 ) + ...    % Middle point proximity
            wt(2) * endError + ...   % Final position accuracy
            wt(3) * timePenalty;     % Time minimization
end



% myTwolinkwithprefilter function
function dxdt= myTwolinkwithprefilter(t, x, t_st, qDes, zeta1,zeta2,zeta3,wn1,wn2,wn3)
    A1 = [zeros(3), eye(3); -diag(wn1).^2, -2 * diag(zeta1) * diag(wn1)];
    B1 = [zeros(3); diag(wn1).^2];

    A2 = [zeros(3), eye(3); -diag(wn2).^2, -2 * diag(zeta2) * diag(wn2)];
    B2 = [zeros(3); diag(wn2).^2];

    A3 = [zeros(3), eye(3); -diag(wn3).^2, -2 * diag(zeta3) * diag(wn3)];
    B3 = [zeros(3); diag(wn3).^2];

    q   = x(7:9);
    qd  = x(10:12);
    
    Kp = diag([140 140 140]);  

    Kd = diag([30 30 30]);  

    controller = Kp * (x(1:3) - q) + Kd * (x(4:6) - qd);
    

    [M, C, G] = compute_M_C_G(q(1), q(2), q(3), qd(1), qd(2), qd(3));
    
    tau = M * (controller) + C * qd ;
    
    qdd = M \ (tau - C * qd );

    if   t <  t_st(1)
        dxdt = [A1*x(1:6) + B1*qDes(:); qd; qdd];
    elseif t >=  t_st(1) && t <  t_st(2)
        dxdt = [A2*x(1:6) + B2*qDes(:); qd; qdd];
    else
        dxdt = [A3*x(1:6) + B3*qDes(:); qd; qdd];
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





