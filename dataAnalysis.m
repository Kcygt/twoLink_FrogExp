% Load data
data = load('The_last_Bump_surface.mat');

% get Position and Force data
step = 2:2500:19001;

pos = [ data.trajectory_first(step,3) data.trajectory_first(step,5)];
fAct = data.trajectory_first(step,12);

desired_size = 9000; % Or 800, or any other desired size
num_knots = size(pos, 1); % Number of points in the input knots
nint = ceil((desired_size - 1) / (num_knots - 1) + 1); % Calculate nint for the desired size

% Use the calculated nint

% Generate the spline with the adjusted size
spline = BSpline(pos, 'order', 4, 'nint', nint);


% spline = BSpline(pos,'order',4,'nint',100);
figure(1); hold on; grid on;
plot(pos(:,1),pos(:,2),'*');
plot(spline(:,1),spline(:,2))

% figure(2); hold on; grid on;
% plot(fAct,'o')


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

	if ip.Results.periodic
		BS=BS(m_min:m_max-1,:);
		BS(end,:)=BS(1,:);
	end
end