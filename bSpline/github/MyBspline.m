close all

% rng(1) % Set random seed for reproductility
% XY=rand(4,2); % Random set of 5 points in 2D

XY = [0	0;
      1 1;
      2 -3;
      3 -2
      ];
BS2=BSpline(XY); % Default order=3 (quadratic)
BS3=BSpline(XY,'order',4); % order=4 -> cubic B-spline
h=plot(XY(:,1),XY(:,2),'-o',BS2(:,1),BS2(:,2),BS3(:,1),BS3(:,2),'--');

legend('Control polygon','Quadratic','Cubic')
set(h, 'LineWidth',2) 


function BS = BSpline(knots,varargin)
%BSPLINE computes the B-spline approximation from a set of coordinates.
% BSPLINE(KNOTS) returns the B-spline interpolation between the reference
% points (knots) whose coordinates are given by the array KNOTS.
% The coordinates of the knots are given vertically, i.e. KNOTS(i,j) gives 
% the j-th coordinate of the i-th knot. The knots can be of any dimension.
%
% BSPLINE(KNOTS,'order',n) uses n -th order approximation, that is 
% polynomial function of degree n-1 (default: n=3, for quadratic B-spline).
%
% BSPLINE(KNOTS,'nint',m) gives m points per interval (default: m=10)
%
% If KNOTS is of size [p,q], the result will be of size [(p-1)*(m-1)+1 ,q],
% except if periodicity is requested (see below).
%
% BSPLINE(KNOTS,'periodic',true) use periodic conditions at end knots.
%

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



% %% Console
% close all;
% clear
% format long;
% clc;
% 
% %% Parameters
% C = [0.0 1.0 4.0 5.0 6 8 10 ;  % x components
%     0.0 1.0 0.5 0.0 -1 -2 1]; % y components
% 
% % B-Spline Degree
% m = 3;
% 
% % B-Spline
% s = f_Bspline(C, m, .1);
% 
% %% Plot
% figure(1);
% plot(C(1,:), C(2,:),'o');
% hold on;
% plot(C(1,:), C(2,:),'--');
% hold on;
% plot(s(1,1:end-1), s(2,1:end-1));
% grid on;
% grid minor;
% xlabel('x');
% ylabel('y');
% legend({'Control points','Polygon','B-Spline-Approx.'}, 'Location', 'southeast');
% ylim([min(C(2,:)) max(C(2,:))*1.25]);
% 
% 
% function s = f_Bspline(C, m, step)
% % Calculates a B-Spline of degree m using the control points C.
% %
% % C: 2-dimensional control points (x_0, ... , x_n; y_0, ... , y_n)
% % m: B-Spline degree
% % step: Step size of t.
% %
% % s: B-Spline. s(1,:) -> x component, s(2,:) -> y component
% 
% %% Parameters
% 
% % Control point's x and y components
% x = C(1,:);
% y = C(2,:);
% 
% % Number of control point - 1
% n = size(x,2) - 1;
% 
% % Knot vector
% T = f_BSpline_KnotVector(m,n);
% 
% % B-Spline intervall
% t = 0:step:(n-m+1);
% 
% %% Calculate B-Spline
% 
% for z=1:1:size(t,2)
%     ti = t(z);
%     s(1,z) = 0;
%     s(2,z) = 0;
%     for i=0:1:n
%         % Base B-Spline
%         B = f_BSpline_BaseSpline(i, m, ti, T);
%         % x component
%         s(1,z) = s(1,z) + x(i+1) * B;
%         % y component
%         s(2,z) = s(2,z) + y(i+1) * B;
%     end
% end
% end
% 
% 
% function T = f_BSpline_KnotVector(m,n)
% % Calculate knot sequence.
% % m: Degree of B-Spline
% % n: Number of control points - 1
% %
% % T = [t0, ... , t_n+m+1] : Knot sequence / vector
% 
% T = zeros(1, (n+m+2));
% for j=0:1:(n+m+1)
%     if j <= m
%         Ti = 0;
%     elseif j >= (m+1) && j <= n
%         Ti = j - m;
%     elseif j > n
%         Ti = n - m + 1;
%     end
%     T(j+1) = Ti;
% end
% 
% end
% 
% 
% function B = f_BSpline_BaseSpline(i, k, t, T)
% % Calculate Base B-Spline B_i,k
% % k: B-Spline degree
% % T: Knot sequence
% % t: Current t parameter
% %
% % B: Base B-Spline at t.
% 
% % Index shift
% j = i + 1;
% 
% if k == 0
%     % End of recusrion
%     if t >= T(j) && t < T(j+1)
%         B = 1;
%     else
%         B = 0;
%     end
% 
% else
%     % Check dividing by zero
%     if T(j+k) ~= T(j)
%         A = (t -T(j))/(T(j+k) - T(j));
%     else
%         A = 0;
%     end
%     if T(j+k+1) ~= T(j+1)
%         B = (T(j+k+1) - t) / (T(j+k+1) - T(j+1));
%     else
%         B = 0;
%     end
%     % Calculate base B-Spline
%     B1 = f_BSpline_BaseSpline(i,   k-1, t, T);
%     B2 = f_BSpline_BaseSpline(i+1, k-1, t, T);
%     B = A * B1 + B * B2;
% end
% end
% 
