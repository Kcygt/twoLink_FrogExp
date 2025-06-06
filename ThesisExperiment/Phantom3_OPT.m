clear; clc;
close all;

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

ttime = [0.8, 1.2, 2.25, 3.75, 5.05 ];
tspan = [0, ttime(end)];

% zeta1 = [.5 1 1]; 
% zeta2 = [0.8 1 0.2];
% zeta3 = [1 1 5];
% zeta4 = [4 1 .7];
% zeta5 = [  1 1 1];
% 
% wn1 = [1 1 1]; 
% wn2 = [1 1 2];
% wn3 = [1 1 1.5];
% wn4 = [1 1  1];
% wn5 = [4 1 4];
zeta1 = [1 1 20]; 
zeta2 = [1 1 10];
zeta3 = [1 1 5];
zeta4 = [1 1 1];
zeta5 = [1 1 1];

wn1 = [1 1 1]; 
wn2 = [1 1 1];
wn3 = [1 1 1];
wn4 = [1 1 1];
wn5 = [1 1 1];

kj1 = [50 50 50];
kj2 = [50 50 50];
kj3 = [50 50 50];
kj4 = [50 50 50];
kj5 = [50 50 50];

bj1 = [30 30 30];
bj2 = [30 30 30];
bj3 = [30 30 30];
bj4 = [30 30 30];
bj5 = [30 30 30];

% weights
wt = [1, .00001, 200];  %  [qDes, Time, qMid]

% Optimization setup
initPrms = [ttime,zeta1,zeta2,zeta3,zeta4,zeta5,wn1,wn2,wn3,wn4,wn5,kj1,kj2,kj3,kj4,kj5,bj1,bj2,bj3,bj4,bj5];

% Initial Condition
[ti, yi] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, initPrms(1:5), qDes, initPrms(6:8),initPrms(9:11),initPrms(12:14),initPrms(15:17),initPrms(18:20), ...   
                                                                            initPrms(21:23),initPrms(24:26),initPrms(27:29),initPrms(30:32),initPrms(33:35), ...
                                                                            initPrms(36:38),initPrms(39:41),initPrms(42:44),initPrms(45:47),initPrms(48:50), ...
                                                                            initPrms(51:53),initPrms(54:56),initPrms(57:59),initPrms(60:62),initPrms(63:65)), ...
                                                                            [0 initPrms(5)], zeros(12, 1));


% Lower and Upper Limits
lb = [0  0  0  0  0   ...                                % time 
      .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .1  ...  % zeta
      .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 ...                  % wn
      70 70 70 70 70 70 70 70 70 70 70 70 70 70 70  ...  % Kp   
      30 30 30 30 30 30 30 30 30 30 30 30 30 30 30  ];   % Kd
ub = [10 10 10 10 10 ...                                     % time
      2  2  2  2  2  2  2  2  5  5  2  2  2  2  2 ...        % zeta
      1  1  1  1  1  2  1  1  2  1  1  1  5  1  5   ...      % wn
      80 80 80 80 80 80 80 80 80 80 80 80 80 80 80 ...       % Kp
      40 40 40 40 40 40 40 40 40 40 40 40 40 40 40  ];       % Kd


% Objective Function
objectiveFunc = @(params) objectiveFunction(params, qDes, wt, qMid,xMid);

% Time contrain Function
constraintFunc = @(prms) timeConstraints(prms);

% Run optimization
options = optimset('PlotFcns','optimplotfval','Display', 'off', 'TolFun', 1e-8, 'MaxIter', 400,'TolX',1e-8);
Opt = fmincon(objectiveFunc, initPrms, [], [], [], [], lb, ub, constraintFunc, options);



% Simulate with optimal parameters
[tt, yy] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, Opt(1:5), qDes, Opt(6:8),   Opt(9:11),  Opt(12:14),  Opt(15:17), Opt(18:20), ...   
                                                                       Opt(21:23), Opt(24:26), Opt(27:29),  Opt(30:32), Opt(33:35), ...
                                                                       Opt(36:38), Opt(39:41), Opt(42:44),  Opt(45:47), Opt(48:50), ...
                                                                       Opt(51:53), Opt(54:56), Opt(57:59),  Opt(60:62), Opt(63:65)), ...
                                                                       [0 Opt(5)], zeros(12, 1));

%%% Plotting
[xi,yi,zi] = FK(yi(:,7),yi(:,8),yi(:,9));  % Initial Trajectory
[x,y,z] = FK(yy(:,7),yy(:,8),yy(:,9));     % Optimized Trajectory

figure; hold on; grid on;
plot(xi,zi,'--')
plot(x,z)
plot(xMid(1,1),xMid(1,3),'*')
plot(xMid(2,1),xMid(2,3),'*')
plot(xMid(3,1),xMid(3,3),'*')
plot(xMid(4,1),xMid(4,3),'*')
plot(0.05,0.05,'o')
legend('Initial Trajectory','Optimized Trajectory')

% figure; hold on; grid on;
% plot(tt, yy(:,1))
% plot(tt, yy(:,3))
% xlabel('Time (s)')
% ylabel('Joint Position (rad)')
% legend('Joint 1','Joint 3')
% % Plot vertical lines
% for t = ttime
%     xline(t, '--k', 'LineWidth', 1.5); % Dashed black line
% end
% 
% figure; hold on;
% colors = lines(length(ttime)); % Generate distinct colors
% for jj = 1:length(ttime)
%     if jj == 1
%         idx = (tt >= 0) & (tt < ttime(jj)); % First stage
%     else
%         idx = (tt >= ttime(jj-1)) & (tt < ttime(jj)); % Subsequent stages
%     end
%     plot(x(idx), z(idx), '-','Color', colors(jj, :), 'LineWidth', 1.5);
% end
% hold off;
% legend(arrayfun(@(jj) sprintf('Stage %d', jj), 1:length(ttime), 'UniformOutput', false));
% xlabel('x');
% ylabel('z');
% title('Plot of different time stages');
% grid on;
% 
% disp('Optimal Parameter:')
disp(['Time: ', num2str(Opt(1:5))])
% disp(['Zeta: ', num2str(Opt(6:20))])
% disp(['Wn: ', num2str(Opt(21:35))])
% disp(['Kp: ', num2str(Opt(36:50))])
% disp(['Kd: ', num2str(Opt(51:65))])



% Constrain Function
function [c, ceq] = timeConstraints(prms)
    % Nonlinear inequality constraint (must be negative or zero)
    c = [prms(1) - prms(2);
         prms(2) - prms(3);
         prms(3) - prms(4);
         prms(4) - prms(5)];
    
    % No equality constraint
    ceq = [];
end



% Objective function
function error = objectiveFunction(prms, qDes, wt, qMid,xMid)
    x0 = zeros(12, 1);
    x0(1:3) = qDes; 

    % Simulate the system
    [t, y] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, prms(1:5), qDes,   prms(6:8),  prms(9:11),  prms(12:14),  prms(15:17), prms(18:20), ...
                                                                              prms(21:23), prms(24:26), prms(27:29),  prms(30:32), prms(33:35), ...
                                                                              prms(36:38), prms(39:41), prms(42:44),  prms(45:47), prms(48:50), ...
                                                                              prms(51:53), prms(54:56), prms(57:59),  prms(60:62), prms(63:65)), ...
                                                                              [0 prms(5)], x0);
    [Cx,~,Cz] = FK(y(:,7),y(:,8),y(:,9));     % Optimized Trajectory
    
    % Calculate error metric
    distto1 = min(sum((y(:, 7:9) - qDes).^2, 2) + sum((prms(5) - t).^2, 2)); 
    distMid1X = min(Cx-xMid(1,1))^2;
    distMid2X = min(Cx-xMid(2,1))^2;
    distMid3X = min(Cx-xMid(3,1))^2;
    distMid4X = min(Cx-xMid(4,1))^2;
    penaltyX = (distMid1X + distMid2X + distMid3X + distMid4X );

    distMid1Z = min(Cz-xMid(1,3) );
    distMid2Z = min(Cz-xMid(2,3) );
    distMid3Z = min(Cz-xMid(3,3) );
    distMid4Z = min(Cz-xMid(4,3) );
    penaltyZ = (distMid1Z + distMid2Z + distMid3Z + distMid4Z );

    % distMid = sum(arrayfun(@(i) min(sum((y(:, 7:9) - qMid(i, :)).^2, 2)), 1:size(qMid,1)));
    error = wt(1)*distto1+wt(3)*(penaltyX + penaltyZ) + wt(2)*prms(5);
    % error = wt(1) * distto1 + wt(2) * prms(5) + wt(3) * distMid;

end


% myTwolinkwithprefilter function
function dxdt= myTwolinkwithprefilter(t, x, t_st, qDes, zeta1,zeta2,zeta3,zeta4,zeta5,wn1,wn2,wn3,wn4,wn5,kj1,kj2,kj3,kj4,kj5,bj1,bj2,bj3,bj4,bj5)


    % Store zeta and wn values in arrays for iteration
    zeta = {zeta1, zeta2, zeta3, zeta4, zeta5};
    wn = {wn1, wn2, wn3, wn4, wn5};
    % kj = {kj1 kj2 kj3 kj4 kj5};
    % bj = {bj1 bj2 bj3 bj4 bj5};

    % Initialize cell arrays for A and B
    A = cell(1, 5);
    B = cell(1, 5);
    
    % Loop through each index and compute A and B matrices
    for i = 1:5
        A{i} = [zeros(3), eye(3); -diag(wn{i}).^2, -2 * diag(zeta{i}) * diag(wn{i})];
        B{i} = [zeros(3); diag(wn{i}).^2];
    end
    q   = x(7:9);
    qd  = x(10:12);
    
    Kp1 = diag(kj1);  
    Kp2 = diag(kj2);  
    Kp3 = diag(kj3);  
    Kp4 = diag(kj4);  
    Kp5 = diag(kj5);  
    
    Kd1 = diag(bj1);  
    Kd2 = diag(bj2);  
    Kd3 = diag(bj3);  
    Kd4 = diag(bj4);  
    Kd5 = diag(bj5);  
    
    controller1 = Kp1 * (x(1:3) - q) + Kd1 * (x(4:6) - qd);
    controller2 = Kp2 * (x(1:3) - q) + Kd2 * (x(4:6) - qd);
    controller3 = Kp3 * (x(1:3) - q) + Kd3 * (x(4:6) - qd);
    controller4 = Kp4 * (x(1:3) - q) + Kd4 * (x(4:6) - qd);
    controller5 = Kp5 * (x(1:3) - q) + Kd5 * (x(4:6) - qd);

    [M, C, G] = compute_M_C_G(q(1), q(2), q(3), qd(1), qd(2), qd(3));
    
    tau1 = M * (controller1) + C * qd ;
    tau2 = M * (controller2) + C * qd ;
    tau3 = M * (controller3) + C * qd ;
    tau4 = M * (controller4) + C * qd ;
    tau5 = M * (controller5) + C * qd ;

    qdd1 = M \ (tau1 - C * qd );
    qdd2 = M \ (tau2 - C * qd );
    qdd3 = M \ (tau3 - C * qd );
    qdd4 = M \ (tau4 - C * qd );
    qdd5 = M \ (tau5 - C * qd );


    if t < t_st(1)
        dxdt = [A{1}*x(1:6) + B{1}*qDes(:); qd; qdd1];
    elseif t_st(1) <= t && t < t_st(2)
        dxdt = [A{2}*x(1:6) + B{2}*qDes(:); qd; qdd2];
    elseif t_st(2) <= t && t < t_st(3)
        dxdt = [A{3}*x(1:6) + B{3}*qDes(:); qd; qdd3];
    elseif t_st(3) <= t && t < t_st(4)
        dxdt = [A{4}*x(1:6) + B{4}*qDes(:); qd; qdd4];
    else
        dxdt = [A{5}*x(1:6) + B{5}*qDes(:); qd; qdd5];
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





