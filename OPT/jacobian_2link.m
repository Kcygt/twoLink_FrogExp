
function J = jacobian_2link(q1, q2, l1, l2)
    J = [
        -l1*sin(q1) - l2*sin(q1 + q2), -l2*sin(q1 + q2);
         l1*cos(q1) + l2*cos(q1 + q2),  l2*cos(q1 + q2)
    ];
end