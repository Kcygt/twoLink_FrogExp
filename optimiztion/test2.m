% Given data
Force_actual = [-0.2264800; 1.0072401; 1.4035800; 1.4721200; 1.5257601; 1.4035800; 0.6615600; 1.0281000];
Force_desired = [1; 1; 1; 1; 1; 1; 1; 1];
position_actual = [0.05 -0.1966063;
                   0.04 -0.1969542;
                   0.03 -0.1905693;
                   0.02 -0.1737738;
                   0.01 -0.1670248;
                   0.00 -0.1694708;
                  -0.01 -0.1856128;
                  -0.02 -0.1981216];

% Parameters for optimization
learning_rate = 0.01; % Step size for position updates
tolerance = 1e-6;     % Convergence tolerance
max_iterations = 100; % Maximum iterations

% Initial z-axis positions
z_current = position_actual(:, 2);

% Iterative optimization
for iter = 1:max_iterations
    % Compute force error
    Force_error = Force_desired - Force_actual;
    
    % Update z-axis positions based on error (proportional control)
    z_next = z_current + learning_rate * Force_error;
    
    % Update actual force using a simple linear model: Force = k * z
    k = 10; % Example proportionality constant (adjust based on your system)
    Force_actual = k * z_next;
    
    % Check for convergence (based on force error)
    if max(abs(Force_error)) < tolerance
        fprintf('Converged after %d iterations.\n', iter);
        break;
    end
    
    % Update current positions
    z_current = z_next;
end

% Display final results
disp('Optimized z-axis positions:');
disp(z_current);
disp('Final Force_actual:');
disp(Force_actual);
figure(1); hold on;
plot(pos(:,1),pos(:,2),'*');
plot(pos(:,1),z_current,'*');
