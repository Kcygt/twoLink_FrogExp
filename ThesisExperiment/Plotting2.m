%%% Plotting
[xi,yi,zi] = FK(yinit(:,7),yinit(:,8),yinit(:,9));  % Initial Trajectory
[x,y,z] = FK(yy(:,7),yy(:,8),yy(:,9));     % Optimized Trajectory

figure; hold on; grid on;
plot(xi,zi,'--')
plot(x,z)
plot(xMid(1,1),xMid(1,3),'*')
plot(xMid(2,1),xMid(2,3),'*')
plot(xMid(3,1),xMid(3,3),'*')
plot(xMid(4,1),xMid(4,3),'*')
plot(0.05,0.0,'o')
legend('Initial Trajectory','Optimized Trajectory')

figure; hold on; grid on;
plot(tt, yy(:,1))
plot(tt, yy(:,3))
xlabel('Time (s)')
ylabel('Joint Position (rad)')
legend('Joint 1','Joint 3')

% Plot vertical lines
for t = ttime
    xline(t, '--k', 'LineWidth', 1.5); % Dashed black line
end

figure; hold on;
colors = lines(length(ttime)); % Generate distinct colors
for jj = 1:length(ttime)
    if jj == 1
        idx = (tt >= 0) & (tt < ttime(jj)); % First stage
    else
        idx = (tt >= ttime(jj-1)) & (tt < ttime(jj)); % Subsequent stages
    end
    plot(x(idx), z(idx), '-','Color', colors(jj, :), 'LineWidth', 1.5);
end
hold off;
legend(arrayfun(@(jj) sprintf('Stage %d', jj), 1:length(ttime), 'UniformOutput', false));
xlabel('x');
ylabel('z');
title('Plot of different time stages');
grid on;

disp('Optimal Parameter:')
disp(['Time: ', num2str(Opt(1:5))])
disp(['Zeta: ', num2str(Opt(6:20))])
disp(['Wn: ', num2str(Opt(21:35))])
% disp(['Kp: ', num2str(Opt(36:50))])
% disp(['Kd: ', num2str(Opt(51:65))])
