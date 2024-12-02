function [q1, q2] = inverse_kinematicsE(x, y, l1, l2)
    % Inverse kinematics with error handling for unreachable targets
    r = sqrt(x^2 + y^2);
    if r > (l1 + l2) || r < abs(l1 - l2)
        warning('Target point (%f, %f) is out of reach', x, y);
        q1 = NaN; q2 = NaN;  % Return NaN for unreachable points
        return;
    end
    cos_q2 = (r^2 - l1^2 - l2^2) / (2 * l1 * l2);
    sin_q2 = sqrt(1 - cos_q2^2); 
    q2 = atan2(sin_q2, cos_q2);
    phi = atan2(y, x);
    psi = atan2(l2 * sin(q2), l1 + l2 * cos(q2));
    q1 = phi - psi;
end