x0 = zeros(8,1);
x0(1:2) = [0.5; 0.4];
[t, y] = ode45(@(t,x) twolinkwithprefilter11(t, x, 10), [0 30], x0);

figure(1); plot(t, y(:, 1:4)); title('Prefilter States');
figure(2); plot(t, y(:, 7:8)); title('Link Dynamics');
aa = forward_kinematics(y(:,5),y(:,6),1,1);
plot(aa(:,1),aa(:,2))

function dxdt = twolinkwithprefilter11(t, x, wn)
    qDes = [ -0.4240,   2.4189;  
          0.1296    1.9552;  
          0.0,      1.5708;  
         -0.5139    1.9552;  
         -0.4240,   2.4189  
        ]; 

    % Prefilter matrix
    zeta = 1; % Adjusted damping ratio
    A = [zeros(2,2), eye(2); -eye(2)*wn^2, -eye(2)*2*zeta*wn];

    % Two-link arm dynamics
    q = x(5:6);
    qp = x(7:8);
    q1p = x(7); q2p = x(8);
    q1 = x(5);  q2 = x(6);

    bj = .01;  % Damping term
    kj = 1.0;  % Control gain
    L_1 = 1; L_2 = 1;
    m_1 = 1; m_2 = 1;

    % Derived constants
    ka = L_2^2 * m_2;
    kb = L_1 * L_2 * m_2;
    kc = L_1^2 * (m_1 + m_2);

    % Mass matrix
    M = [ka + 2*kb*cos(q2) + kc, ka + kb*cos(q2);
         ka + kb*cos(q2), ka];
    C = coriolis_matrix(q1, q2, q1p, q2p, L_1, L_2, m_1, m_2);
    
    % Control law
    Numerator = [-bj, 0; 0, -bj] * qp  + [-kj, 0; 0, -kj] * (q - x(1:2)) - C*qp;
    qdd = M \ Numerator;

    % Prefilter dynamics
    if t < 0.4
        dotx = [0; 0; 0; 0];
    elseif t < 5
        Uinp = [.2; 0.3; 0; 0] * wn^2;
        dotx = A * x(1:4) + Uinp;
    else
        Uinp = [0; 0; 0.0; 0] * wn^2;
        dotx = A * x(1:4) + Uinp;
    end

    dxdt = [dotx; qp; qdd];
end
