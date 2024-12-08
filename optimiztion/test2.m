% Load data
data = load('The_last_Bump_surface.mat');

% Get Position and Force data
step = 2:2500:19001;

pos = [data.trajectory_first(step,3), data.trajectory_first(step,5)];
fAct = data.trajectory_first(step,12);

desired_size = 1000; % Desired total number of points
num_knots = size(pos, 1); % Number of input knots

% Uniform parameterization
nint = ceil(desired_size / (num_knots - 1)); % Points per segment
total_points = (num_knots - 1) * nint + 1; % Recalculate total points
spline = UniformBSpline(pos, 'order', 4, 'nint', nint);

% Plot results
figure(1); hold on; grid on;
plot(pos(:,1), pos(:,2), '*'); % Original knots
plot(spline(1:125,1), spline(1:125,2), '-'); % Interpolated spline



function BS = UniformBSpline(knots, varargin)
    ip = inputParser;
    addOptional(ip, 'order', 3); % B-spline order
    addOptional(ip, 'nint', 10); % Points per interval
    addOptional(ip, 'periodic', false); % Periodicity
    parse(ip, varargin{:});
    
    % Handle periodicity
    if ip.Results.periodic
        np_rep = ip.Results.order;
        knots = [knots(end-np_rep+1:end, :); knots; knots(1:np_rep, :)];
    end
    
    % Dimensions
    p = size(knots, 1); % Number of input points
    q = size(knots, 2); % Dimensionality of data
    
    if p <= 2
        BS = knots; % Too few points, return input
        return;
    end
    
    nint = ip.Results.nint; % Points per segment
    t = linspace(0, 1, nint * (p - 1) + 1); % Uniform parameterization
    order = min(ip.Results.order, p); % Spline order
    
    % Knot vector for clamped spline
    tk = [zeros(1, order-1), linspace(0, 1, p-order+2), ones(1, order-1)];
    
    % Initialize output
    BS = zeros(length(t), q);
    basis = zeros(order, q); % Initialize basis array for dimensions

    % De Boor's algorithm
    for i = 1:length(t)
        t0 = t(i);
        k = find(t0 >= tk, 1, 'last');
        k = min(k, p); % Clamp to valid range

        % Initialize basis with the correct dimensionality
        for d = 1:q % Loop over dimensions
            basis(:, d) = knots(max(1, k-order+1):min(p, k), d);
        end

        for j = 2:order
            for l = j:order
                alpha = (t0 - tk(k-order+l)) / (tk(k+l-j+1) - tk(k-order+l));
                alpha(isnan(alpha) | isinf(alpha)) = 0; % Handle edge cases
                for d = 1:q % Loop over dimensions
                    basis(l, d) = (1-alpha) * basis(l-1, d) + alpha * basis(l, d);
                end
            end
        end
        BS(i, :) = basis(order, :);
    end
    
    % Handle periodicity output
    if ip.Results.periodic
        BS = BS(1:end-nint+1, :);
        BS(end, :) = BS(1, :);
    end
end
