function C = coriolis_matrix(q1, q2, q1d, q2d, l1, l2, m1, m2)
    C1 = -l1*l2*m2*q2d*sin(q2)*(2*q1d + q2d);
    C2 = l1*l2*m2*q1d^2*sin(q2);
    C = [C1; C2];
end