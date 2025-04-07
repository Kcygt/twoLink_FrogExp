% Output for initial contiditions
[xInit, yInit, zInit] = FK(init_Y(:,7), init_Y(:,8), init_Y(:,9));



% Desired Trajectory for initial condition
[xD, yD, zD] = FK(y(:,1), y(:,2), y(:,3));

% Output
[xAct, yAct, zAct] = FK(y(:,7), y(:,8), y(:,9));
[xDes, yDes, zDes] = FK(qDes(1), qDes(2), qDes(3));

% Plotting
figure; hold on; grid on;
plot(xInit, zInit, '-.');
plot(xAct, zAct, '.-');
plot(xD,zD);
plot(xDes, zDes, 'o');
plot(xMid(:,1),xMid(:,3),'*')
xlabel('X axis'); ylabel('Z axis');
legend('Initial', 'Optimized', 'Desired trajectory','Final Points','Mid points')
title('Cartesian Trajectory Tracking');
% 
% figure; hold on; grid on;
% plot(y(:,7),y(:,9),'-')
% plot(qMid(:,1),qMid(:,3),'o')
% xlabel('Joint 1')
% ylabel('Joint 2')
% title('Joint Space Trajectory')
% 
% 
% figure; hold on; grid on;
% plot(y(:,7))
% plot(y(:,1))
% plot(qMid(:,1),'*');
% legend('Actual Joint1 position','Desired Joint1 Position','Mid-Points')
% 
% 
disp(['Optimized Parameters: ', num2str(optimalParams)]);
