function M = mass_matrix(q1, q2, l1, l2, m1, m2)
    M11 = l1^2*m1 + l1^2*m2 + l2^2*m2 + 2*l1*l2*m2*cos(q2);
    M12 = l2*m2*(l2 + l1*cos(q2));
    M22 = l2^2*m2;
    M = [M11, M12; M12, M22];
end