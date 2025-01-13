function G = gravity_vector(q1, q2, l1, l2, m1, m2, gx, gy)
    G1 = gx * (l1 * m1 * sin(q1) + l1 * m2 * sin(q1) + l2 * m2 * sin(q1 + q2)) + gy * (l1 * m1 * cos(q1) + l1 * m2 * cos(q1) + l2 * m2 * cos(q1 + q2));
    G2 = gx * l2 * m2 * sin(q1 + q2) + gy * l2 * m2 * cos(q1 + q2);
    G = [G1; G2];
end
