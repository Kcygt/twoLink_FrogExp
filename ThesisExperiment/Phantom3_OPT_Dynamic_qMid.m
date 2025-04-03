clear; clc;
close all;

% Define desired trajectory and Middle Points
qDes = [0.114, -0.0445, 0.3336];
xMid = [0.03, 0, 0.01 ];
qMid = IK(xMid(1), xMid(2), xMid(3));

% Parameters
ttime = [5 10];
tspan = [0, ttime(end)];

zeta1 = [.5 1 0.1]; 
zeta2 = [0.8 1 0.1];

wn1 = [10 1 1]; 
wn2 = [1 1 20];

kj1 = [50 50 50];
kj2 = [50 50 50];

bj1 = [30 30 30];
bj2 = [30 30 30];

% weights
% wt = [1, 1e-10, 30,100,1e+3];  %  [qDes, Time, qMid,xMid,Velocity]
wt = [0 0 0 0 1e+3];  %  [qDes, Time, qMid,xMid,Velocity]

% Optimization setup
initPrms = [ttime,zeta1,zeta2,wn1,wn2,kj1,kj2,bj1,bj2,qMid];

% Initial Condition
[tInitial, yInitial] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, initPrms(1:2), qDes, initPrms(3:5),   initPrms(6:8),   ...
                                                                            initPrms(9:11),  initPrms(12:14), ...
                                                                            initPrms(15:17), initPrms(18:20), ...   
                                                                            initPrms(21:23), initPrms(24:26), initPrms(27:29)), ...
                                                                            [0 initPrms(2)], zeros(12, 1));

% Lower and Upper Limits
lb = [ 0  0     ...          % time 
      .1 .1 .1 .1 .1 .1 ...  % zeta
      .5 .5 .5 .5 .5 .5 ...  % wn
      70 70 70 70 70 70 ...  % Kp   
      30 30 30 30 30 30 ...
      0 0 0];   % Kd

ub = [10 10  ...              % time
      1  1  1  1  1  1  ...   % zeta
      20  20  20  20  20  20  ...   % wn
      80 80 80 80 80 80 ...   % Kp
      40 40 40 40 40 40 ... 
      0.5 0.5 0.5];    % Kd


% Objective Function
objectiveFunc = @(params) objectiveFunction(params, qDes, wt);

% Time contrain Function
constraintFunc = @(prms) timeConstraints(prms);

% Run optimization
options = optimset('PlotFcns','optimplotfval','Display', 'off', 'TolFun', 1e-8, 'MaxIter', 400,'TolX',1e-8);
OptParams = fmincon(objectiveFunc, initPrms, [], [], [], [], lb, ub, constraintFunc, options);


% Simulate with optimal parameters
[tOutput, yOutput] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, OptParams(1:2), qDes, OptParams(3:5),   OptParams(6:8),   ...
                                                                             OptParams(9:11),  OptParams(12:14), ...
                                                                             OptParams(15:17), OptParams(18:20), ...
                                                                             OptParams(21:23), OptParams(24:26),...
                                                                             OptParams(27:29)), ...
                                                                             [0 OptParams(2)], zeros(12, 1));

%%% Plotting
[xi,~,zi] = FK(yInitial(:,7),yInitial(:,8),yInitial(:,9));  % Initial Trajectory
[xOpt,yOpt,zOpt] = FK(yOutput(:,7),yOutput(:,8),yOutput(:,9));     % Optimized Trajectory

[xMidO,yMidO,zMidO] = FK(OptParams(27),OptParams(28),OptParams(29));     % Optimized Trajectory

figure; hold on; grid on;
plot(xi,zi,'--')     % Initial condition
plot(xOpt,zOpt)            % Optimal condition
plot(xMidO,zMidO,'*')
plot(0.05,0.05,'o')
legend('Initial Trajectory','Optimized Trajectory')

disp('Optimal Parameter:')
disp(['Time: ', num2str(OptParams(1:2))])
disp(['Zeta: ', num2str(OptParams(3:8))])
disp(['Wn: ', num2str(OptParams(9:14))])
disp(['Kp: ', num2str(OptParams(15:20))])
disp(['Kd: ', num2str(OptParams(21:26))])
disp(['Velocity: ', num2str(yOutput(end,10:12))])




% Constrain Function
function [c, ceq] = timeConstraints(prms)
    % Nonlinear inequality constraint (must be negative or zero)
    c = prms(1) - prms(2);
    
    % No equality constraint
    ceq = [];
end



% Objective function
function error = objectiveFunction(prms, qDes, wt)
    x0 = zeros(12, 1);
    x0(1:3) = qDes; 

    % Simulate the system
    [t, y] = ode23s(@(t, x) myTwolinkwithprefilter(t, x, prms(1:2), qDes,   prms(3:5),   prms(6:8), ...
                                                                            prms(9:11),  prms(12:14), ...
                                                                            prms(15:17), prms(18:20), ...
                                                                            prms(21:23), prms(24:26),prms(27:29)), ...
                                                                            [0 prms(2)], x0);
    
    
    xMid = [ 0.01,0,0.01; 
             0.02,0,0.02;
             0.03,0,0.03;
             0.04,0,0.04];
    
    tMid = linspace(0,max(t),6);
    yyintp=interp1(t,y,tMid(2:end-1));
    [xAct,yAct,zAct] = FK(yyintp(:,7),yyintp(:,8),yyintp(:,9));
    

  
    distErr1 = sqrt( (xMid(1,1) - xAct(1)).^2  +  (xMid(1,3) - zAct(1)).^2 ) ;
    distErr2 = sqrt( (xMid(2,1) - xAct(2)).^2  +  (xMid(2,3) - zAct(2)).^2 ) ;
    distErr3 = sqrt( (xMid(3,1) - xAct(3)).^2  +  (xMid(3,3) - zAct(3)).^2 ) ;
    distErr4 = sqrt( (xMid(4,1) - xAct(4)).^2  +  (xMid(4,3) - zAct(4)).^2 ) ;

    % Calculate error metric
    distto1 = min(sum((y(:, 7:9) - qDes).^2, 2) + sum((prms(2) - t).^2, 2)); 

    distMid = sum(arrayfun(@(i) min(sum((y(:, 7:9) - prms(27:29)).^2, 2)), 1:1));
    xError = distErr1 + distErr2 + distErr3 + distErr4;
    Vel = abs(y(end,10)+y(end,11)+y(end,12));
    error = wt(1) * distto1 + wt(2) * prms(2) + wt(3) * distMid + wt(4)*(xError) + wt(5) * Vel;


end


% myTwolinkwithprefilter function
function dxdt= myTwolinkwithprefilter(t, x, t_st, qDes, zeta1,zeta2,wn1,wn2,kj1,kj2,bj1,bj2,qMid)


    % Store zeta and wn values in arrays for iteration
    zeta = {zeta1, zeta2};
    wn =   {wn1, wn2};
    kj = {kj1 kj2};
    bj = {bj1 bj2};

    % Initialize cell arrays for A and B
    A = cell(1, 2);
    B = cell(1, 2);
    
    % Loop through each index and compute A and B matrices
    for i = 1:2
        A{i} = [zeros(3), eye(3); -diag(wn{i}).^2, -2 * diag(zeta{i}) * diag(wn{i})];
        B{i} = [zeros(3); diag(wn{i}).^2];
    end
    q   = x(7:9);
    qd  = x(10:12);
    
    Kp1 = diag(kj{1});  
    Kp2 = diag(kj{2});  
  
    Kd1 = diag(bj{1});  
    Kd2 = diag(bj{2});  
  
    controller1 = Kp1 * (x(1:3) - q) + Kd1 * (x(4:6) - qd);
    controller2 = Kp2 * (x(1:3) - q) + Kd2 * (x(4:6) - qd);
   
    [M, C, G] = compute_M_C_G(q(1), q(2), q(3), qd(1), qd(2), qd(3));
    
    tau1 = M * (controller1) + C * qd ;
    tau2 = M * (controller2) + C * qd ;
    
    qdd1 = M \ (tau1 - C * qd );
    qdd2 = M \ (tau2 - C * qd );
    
    if t < t_st(1)
        dxdt = [A{1}*x(1:6) + B{1}*qDes(:); qd; qdd1];
    else
        dxdt = [A{2}*x(1:6) + B{2}*qDes(:); qd; qdd2];
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





