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

% Initial guess for optimization (use the current z-axis values)
z_initial = position_actual(:, 2);

% Define the cost function for optimization
cost_function = @(z) sum((Force_desired - computeForce(z)).^2);

% Define bounds (if needed, e.g., physical constraints)
lb = -0.5 * ones(size(z_initial));  % Example lower bounds
ub = 0.5 * ones(size(z_initial));   % Example upper bounds

% Set up options for fmincon
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% Run optimization
z_optimized = fmincon(cost_function, z_initial, [], [], [], [], lb, ub, [], options);

% Display the results
disp('Optimized z-axis positions:');
disp(z_optimized);
figure(1); hold on;
plot(pos(:,1),pos(:,2),'*');
plot(pos(:,1),z_optimized,'*');


% Function to compute force based on position
function Force = computeForce(z)
    % Example model: Force = k * z; Adjust 'k' based on your system
    k = 100;  % Adjust this constant as per your system's calibration
    Force = k * z;
end
