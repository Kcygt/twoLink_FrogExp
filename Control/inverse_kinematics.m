function qDes = inverse_kinematics(x, y, l1, l2)
    r = sqrt(x^2 + y^2);
    if r > (l1 + l2) || r < abs(l1 - l2)
        error('Target point is out of reach');
    end
    cos_q2 = (r^2 - l1^2 - l2^2) / (2 * l1 * l2);
    sin_q2 = sqrt(1 - cos_q2^2); 
    q2 = atan2(sin_q2, cos_q2);
    phi = atan2(y, x);
    psi = atan2(l2 * sin(q2), l1 + l2 * cos(q2));
    q1 = phi - psi;
    qDes = [q1, q2];
end