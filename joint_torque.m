function tau = joint_torque(q, qd, qdd, l1, l2, m1, m2, g)
    q1 = q(1); q2 = q(2);
    q1d = qd(1); q2d = qd(2);
    q1dd = qdd(1); q2dd = qdd(2);
    M = mass_matrix(q1, q2, l1, l2, m1, m2);
    C = coriolis_matrix(q1, q2, q1d, q2d, l1, l2, m1, m2);
    G = gravity_vector(q1, q2, l1, l2, m1, m2, g);
    tau = M * [q1dd; q2dd] + C + G;
end
