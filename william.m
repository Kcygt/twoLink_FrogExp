% odefun=@(t,x) mysf(t,x,l1,l2,m1,m2,K,B)
% [t,y] = ode45(odefun,[0 20], [0;0;0;0])

% use interp1 to get equal spaced times

function dx = mysf(t,x,l1,l2,m1,m2,K,B)
%MYSF Summary of this function goes here
%   Detailed explanation goes here

    qAct=x(1:2);
    qdAct=x(3:4);

    % Compute dynamics
    M = mass_matrix(qAct(1), qAct(2), l1, l2, m1, m2);
    G = gravity_vector(qAct(1), qAct(2), l1, l2, m1, m2, 0, 0);
    C = coriolis_matrix(qAct(1), qAct(2), qdAct(1), qdAct(2), l1, l2, m1, m2);
    
    % PD control with gravity compensation

    % Error in position and velocity
     
    e =  qAct - [1;1];
    eDot = qdAct - [0;0] ;


    % Torque = K * -e * exp(-e'*0*e) + B * -eDot;
    Torque = K * -e  + B * -eDot;
  
    qddAct = M \ (Torque - C(:)) ;% Update velocity
    % return dx
    dx=[qdAct(:);qddAct(:)];

end

