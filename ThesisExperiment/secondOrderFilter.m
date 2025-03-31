clc; clear; 
% close all;

% Define system parameters
zeta =4;      % Damping ratio
wn = 1;          % Natural frequency (rad/s)

% Second-order system state-space representation
A2 = [0 1; -wn^2 -2*zeta*wn];
B2 = [0; wn^2];
C2 = [1 0];  % Output is displacement (x)
D2 = 0;

% Create state-space system
sys2 = ss(A2, B2, C2, D2);

% Plot step response
figure(3);
step(sys2);
xlim([0 1.75])
title('Step Response of Second-Order System');
xlabel('Time (s)');
ylabel('Response');
grid on;


% % Third-order system state-space representation
% A3 = [0 1 0; 0 0 1; -wn^3 -3*zeta*wn^2 -3*zeta^2*wn];
% B3 = [0; 0; wn^3];
% C3 = [1 0 0]; % Output is displacement (x)
% D3 = 0;
% 
% 
% % Plot Bode plot
% figure;
% bode(sys2);
% title('Bode Plot of Second-Order System');
% grid on;

% 
% % Create state-space system
% sys3 = ss(A3, B3, C3, D3);
% 
% % Plot step response
% figure;
% step(sys3);
% title('Step Response of Third-Order System');
% xlabel('Time (s)');
% ylabel('Response');
% grid on;
% 
% 
% 
% % Plot Bode plot
% figure;
% bode(sys3);
% title('Bode Plot of Third-Order System');
% grid on;