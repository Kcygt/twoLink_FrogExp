function dxdt = robot_dynamics(t, x, l1, l2, m1, m2, g, K, B)
    % Unpack state variables
    q = x(1:2);
    qd = x(3:4);

    % Desired trajectory
    q_des = [pi/4; pi/4]; % Fix for simplicity
    qd_des = [0; 0];

    % Error in position and velocity
    e = q - q_des;
    eDot = qd - qd_des;

    % Dynamics
    M = mass_matrix(q(1), q(2), l1, l2, m1, m2);
    G = gravity_vector(q(1), q(2), l1, l2, m1, m2, 0,0);
    C = coriolis_matrix(q(1), q(2), qd(1), qd(2), l1, l2, m1, m2);

    % PD control with gravity compensation
    Torque = -K * e - B * eDot;

    % State derivatives
    qdd = M \ (Torque - C * qd - G); % Joint accelerations
    dxdt = [qd; qdd];
end
