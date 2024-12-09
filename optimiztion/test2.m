% Load data
data = load('The_last_Bump_surface.mat');

% Get Position and Force data
step = 2:2500:19001;
pos = [data.trajectory_first(step, 3), data.trajectory_first(step, 5)];
fAct = data.trajectory_first(step, 12);
fDes = 1;

% Ensure data is of type double
pos = double(pos);
fAct = double(fAct);
fDes = double(fDes);

% Optimization setup
w = 0.005; % Weight parameter
lb = double(pos(:, 2) - 0.1); % Lower bound for position adjustment
ub = double(pos(:, 2) + 0.1); % Upper bound for position adjustment
opt_pos = pos; % Initialize optimized positions

% Objective function for fmincon
objective = @(adjustment) sum((fDes - (fAct + w * adjustment)).^2);

% Optimize each position point independently
for i = 1:length(pos)
    % Initial guess for optimization
    adjustment0 = 0;

    % Bounds for adjustment
    lb_i = double(lb(i) - pos(i, 2));
    ub_i = double(ub(i) - pos(i, 2));

    % Gradient influence (optional)
    if i > 1
        grad_f = fAct(i) - fAct(i - 1);
    else
        grad_f = 0;
    end

    % Dynamic weight
    dynamic_w = abs(w * (fDes - abs(fAct(i))^2));

    % Perform optimization using fmincon
    adjustment = fmincon(@(adj) objective(adj) + 0.5 * dynamic_w * abs(grad_f) * adj^2, ...
                         adjustment0, [], [], [], [], lb_i, ub_i);

    % Update optimized position
    opt_pos(i, 2) = pos(i, 2) + adjustment;
end

% Generate spline for visualization
desired_size = 1000;
num_knots = size(pos, 1);
nint = ceil((desired_size - 1) / (num_knots - 1) + 1);

spline = BSpline(pos, 'order', 3, 'nint', nint);
spline_opt = BSpline(opt_pos, 'order', 3, 'nint', nint);

% Plot results
figure(1); hold on; grid on;
plot(pos(:, 1), pos(:, 2), '*', 'DisplayName', 'Original Points');
plot(spline(:, 1), spline(:, 2), 'DisplayName', 'Original Spline');
plot(opt_pos(:, 1), opt_pos(:, 2), 'o', 'DisplayName', 'Optimized Points');
plot(spline_opt(:, 1), spline_opt(:, 2), 'DisplayName', 'Optimized Spline');
legend;

% BSpline function remains unchanged

function BS = BSpline(knots,varargin)
	ip = inputParser;
	addOptional(ip,'order',3)
	addOptional(ip,'nint',10)
	addOptional(ip,'periodic',false)
	parse(ip,varargin{:});
	
	if ip.Results.periodic
		np_rep=ip.Results.order;
		knots=[knots(end-np_rep+1:end,:); knots; knots(1:np_rep,:)];
	end	
	
	p=size(knots,1);
	q=size(knots,2);
	
 	if p<=2
 		BS=knots;
 		return
	end
	
	n=ip.Results.nint;
	n=(n-1)*(p-1)+1;	% Overall number of queried points
	y = linspace(0,1,n);

	order=min(ip.Results.order,p);
	Xl = zeros(order,order,q);
	t = [zeros(1,order-1),linspace(0,1,p-order+2),ones(1,order-1)];

	BS=zeros(n,q);
	m_min=1;
	m_max=n;
	for m = 1:n-1
		t0 = y(m);
		k = find(t0 >= t,1,'last');
		if (k > p)
			BS=BS(1:m-1,:);
			return;
		end
		Xl(:,1,:) = knots(k-order+1:k,:);
		if k<=order+1
			m_min=max(m,m_min);
		end
		if k>=p-order+2
			m_max=min(m,m_max);
		end

		for i = 2:order
			for j = i:order
				num = t0-t(k-order+j);
				if num == 0
					wt = 0;
				else
					s = t(k+j-i+1)-t(k-order+j);
					wt = num/s;
				end
				Xl(j,i,:) = (1-wt)*Xl(j-1,i-1,:) + wt*Xl(j,i-1,:);
			end
		end
		BS(m,:) = Xl(order,order,:);
	end
	BS(end,:)=knots(end,:);

	% if ip.Results.periodic
	% 	BS=BS(m_min:m_max-1,:);
	% 	BS(end,:)=BS(1,:);
	% end
end