% Data initialization
time = 0:9; % Time steps
surface_x = [0 1 2 3 4 5 6 7 8 9]; % Surface x-coordinates
surface_y = [0 1 -1 2 3 1 0 2 3 0]; % Surface y-coordinates
force = [1 1 0 3 5 5 1 1 2 0]; % Measured forces along y
path_x = [0 1 2 3 4 5 6 7 8 9]; % Path x-coordinates
path_y = [0 1 2 1 -1 -2 0 2 1 1]; % Initial path y-coordinates

desired_force = 1; % Desired force in N
optimized_path_y = path_y; % Initialize the optimized path

% Optimization settings
for t = 2:length(time)-2
    % Current, future, and temporary data
    curr_y = path_y(t);        % Current y-coordinate
    future_y = path_y(t+2);    % y-coordinate at t+2
    temp_y = path_y(t+1);      % Temporary y-coordinate at t+1 (initial guess)
    curr_f = force(t);         % Current force
    
    % Nested optimization for time t+1
    objective = @(temp_y_new) abs(desired_force - ...
        estimate_force(curr_y, temp_y_new, future_y, curr_f));
    
    % Constraints for the temporary y-coordinate
    lb_temp = temp_y - 2; % Lower bound for temporary y
    ub_temp = temp_y + 2; % Upper bound for temporary y
    temp_y_new = fmincon(objective, temp_y, [], [], [], [], lb_temp, ub_temp);
    
    % Final optimization for the intermediate point
    objective_final = @(y_opt) abs(desired_force - ...
        estimate_force(curr_y, y_opt, future_y, curr_f));
    
    lb = curr_y - 2; % Lower bound for y at t+1
    ub = curr_y + 2; % Upper bound for y at t+1
    y_optimized = fmincon(objective_final, temp_y_new, [], [], [], [], lb, ub);
    
    % Update the optimized path
    optimized_path_y(t+1) = y_optimized;
end

% Display the optimized path
disp('Optimized Path along Y:');
disp(optimized_path_y);

% Plotting
figure;
hold on;

% Plot the surface
plot(surface_x, surface_y, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Surface');

% Plot the non-optimized path
plot(path_x, path_y, 'r-o', 'LineWidth', 1.5, 'DisplayName', 'Non-Optimized Path');

% Plot the optimized path
plot(path_x, optimized_path_y, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Optimized Path');

% Add labels and legend
xlabel('X-Coordinate');
ylabel('Y-Coordinate');
title('Comparison of Paths and Surface');
legend('Location', 'best');
grid on;
hold off;

% Helper function to estimate force based on path
function f_estimated = estimate_force(y0, y1, y2, f0)
    % Example force model using interpolation
    % Replace this with your actual force-surface interaction model
    delta_y1 = y1 - y0;
    delta_y2 = y2 - y1;
    f_estimated = f0 + 0.5 * (delta_y1 + delta_y2); % Simple linear model
end



% % Data initialization (same as above)
% time = 0:10; % Time steps
% surface_x = [0 1 2 3 4 5 6 7 8 9 10]; % Surface x-coordinates
% surface_y = [0 1 1 2 2 3 2 2 1 1 0]; % Surface y-coordinates
% force     = [1 1 0 3 5 5 1 1 2 0 1]; % Measured forces along y
% path_x    = [0 1 2 3 4 5 6 7 8 9 10]; % Path x-coordinates
% path_y    = [0 1 2 1 -1 -2 0 2 1 1 1]; % Initial path y-coordinates
% 
% desired_force = 1; % Desired force in N
% optimized_path_y = path_y; % Initialize the optimized path
% 
% % Optimization settings
% for t = 2:length(time)-2
%     % Current and next two points (time t, t+1, t+2)
%     curr_t = time(t);       % Current time
%     next_t1 = time(t+1);    % Next time step
%     next_t2 = time(t+2);    % Next time step after that
% 
%     % Force and path data for optimization
%     curr_y = path_y(t);     % Current y-coordinate
%     next_y1 = path_y(t+1);  % Next y-coordinate
%     next_y2 = path_y(t+2);  % y-coordinate after next
% 
%     curr_f = force(t);      % Current force
%     next_f1 = force(t+1);   % Next force
%     next_f2 = force(t+2);   % Force after next
% 
%     % Optimization function to minimize force error
%     objective = @(y_new) abs(desired_force - interpolate_force(curr_t, next_t1, next_t2, curr_y, next_y1, y_new));
% 
%     % Constraints (if any)
%     lb = curr_y - 2; % Lower bound for y (example constraint)
%     ub = curr_y + 2; % Upper bound for y (example constraint)
% 
%     % Solve for optimized point
%     y_opt = fmincon(objective, next_y1, [], [], [], [], lb, ub);
% 
%     % Update the path
%     optimized_path_y(t+1) = y_opt;
% end
% 
% % Display the optimized path
% disp('Optimized Path along Y:');
% disp(optimized_path_y);
% 
% % Plotting
% figure;
% hold on;
% 
% % Plot the surface
% plot(surface_x, surface_y, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Surface');
% 
% % Plot the non-optimized path
% plot(path_x, path_y, 'r-o', 'LineWidth', 1.5, 'DisplayName', 'Non-Optimized Path');
% 
% % Plot the optimized path
% plot(path_x, optimized_path_y, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Optimized Path');
% 
% % Add labels and legend
% xlabel('X-Coordinate');
% ylabel('Y-Coordinate');
% title('Comparison of Paths and Surface');
% legend('Location', 'best');
% grid on;
% hold off;
% 
% % Helper function to estimate force based on path
% function f_estimated = interpolate_force(t0, t1, t2, y0, y1, y2)
%     % Example force model using linear interpolation
%     % Replace this with your actual force-surface interaction model
%     f_estimated = (y1 + y2 - y0) / (t1 - t0 + t2 - t1); 
% end
