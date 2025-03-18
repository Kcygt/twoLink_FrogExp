clc; clear; 

% Define parameters for each joint (3-DOF system)
zeta = [1 1 1];  % Damping ratios for each joint
wn   = [1 1 1];  % Natural frequencies
a    = 1;  % Gain factor (can be adjusted)

% Define A, B, C, D for the state-space system
A1 = [  0   1   0   0   0   0   0   0   0;
       0   0   1   0   0   0   0   0   0;
    -a*wn(1)^2  -(a+2*zeta(1)*wn(1))  -(wn(1)^2 + 2*zeta(1)*wn(1)^2)  0   0   0   0   0   0;
       0   0   0   0   1   0   0   0   0;
       0   0   0   0   0   1   0   0   0;
       0   0   0   -a*wn(2)^2  -(a+2*zeta(2)*wn(2))  -(wn(2)^2 + 2*zeta(2)*wn(2)^2)  0   0   0;
       0   0   0   0   0   0   0   1   0;
       0   0   0   0   0   0   0   0   1;
       0   0   0   0   0   0   -a*wn(3)^2  -(a+2*zeta(3)*wn(3))  -(wn(3)^2 + 2*zeta(3)*wn(3)^2) ];

B1 = [ 0  0  0;
      0  0  0;
      a*wn(1)^2  0  0;
      0  0  0;
      0  0  0;
      0  a*wn(2)^2  0;
      0  0  0;
      0  0  0;
      0  0  a*wn(3)^2 ];

C1 = [1 0 0  0 0 0  0 0 0;  % Output for joint 1
     0 0 0  1 0 0  0 0 0;  % Output for joint 2
     0 0 0  0 0 0  1 0 0]; % Output for joint 3

D1 = zeros(3,3);

% Create state-space system
sys = ss(A1, B1, C1, D1);

% Step response for each joint
figure;
step(sys)
title('Step Response of the Third-Order Prefilter (3-DOF System)')
grid on;


% zeta = .8;
% wn   = .5;
% a    = 10;
% A = [ 0 1 0; 0 0 1; -a*wn^2 -(a+2*zeta*wn) -(wn^2 + 2*zeta*wn^2) ];
% B = [ 0; 0; a*wn^2];
% C = [1 0 0];
% D = 0;
% sys = ss(A, B, C, D);  % State-space to transfer function
% step(sys)


% wn = [5 6 7];  % Natural frequency
% zeta = [ 2 3 4];  % Damping ratios
% 
% 
% A = [zeros(3,3), eye(3);
%      -diag(wn).^2, -2*diag(zeta)*diag(wn)];
% B = [zeros(3,3); diag(wn).^2];
% C = eye(3);
% D = zeros(3,3);
% 
% 
% sys = ss(A, B, C, D);  % State-space to transfer function


% % Given parameters
% wn = 2;  % Natural frequency
% zeta = 1;  % Damping ratios
% 
% % Define the diagonal transfer function matrix
% G = tf(9, [1, 2*zeta*wn, wn^2]);


% 
% % Define the natural frequency (omega_n) and damping ratio (zeta)
% wn = 5;     % Example natural frequency
% zeta = 1;    % Example damping ratio
% % K = 4;         % Example gain
% 
% % Define the A, B, C, and D matrices for the state-space representation
% A = [0, 1, 0;
%      0, 0, 1;
%      -wn^2, -2*zeta*wn, -wn^2];
% B = [0; 0; wn^2];
% C = [1, 0, 0];
% D = 0;
% 
% % Create the state-space system
% sys = ss(A, B, C, D);
% 
% step(sys); hold on




%%%%%%%%%%%%%   other form
% Gain = 2;
% % Define the coefficients of the third-order system
% a2 = 5;  % Example coefficient for s^2 term
% a1 = 2;  % Example coefficient for s term
% a0 = Gain;  % Example coefficient for constant term
% K = Gain;   % Example gain
% 
% % Define the A, B, C, and D matrices for the state-space representation
% A = [0, 1, 0;
%      0, 0, 1;
%      -a0, -a1, -a2];
% B = [0; 0; K];
% C = [1, 0, 0];
% D = 0;
% 
% % Create the state-space system
% sys = ss(A, B, C, D);
% step(sys); hold on






% %%%%%%%
% % Define the natural frequency (omega_n) and damping ratio (zeta)
% omega_n = 10;  % Example natural frequency
% zeta = 1;   % Example damping ratio
% 
% % Define the A, B, C, and D matrices for the state-space representation
% A = [0, 1, 0;
%      0, 0, 1;
%      -omega_n^2, -2*zeta*omega_n, -omega_n^2];
% B = [0; 0; omega_n^2];
% C = [1, 0, 0];
% D = 0;
% 
% % Create the state-space system
% sys = ss(A, B, C, D);
% 
% step(sys)