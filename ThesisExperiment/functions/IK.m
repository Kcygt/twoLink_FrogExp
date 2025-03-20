function [q1,q2,q3] = IK(x,y,z)
    l1 = 0.208; 
    l2 = 0.168;  
    q1 = atan2(x,z+l1);
    
    R = sqrt(x^2 + (z+l1)^2);
    r = sqrt(x^2 + (y-l2)^2 + (z+l1)^2);
    Beta  = atan2(y-l2,R);
    Gamma =  acos((l1^2+r^2 - l2^2)/(2*l1*r));
    q2 = Gamma + Beta;

    Alpha = acos((l1^2 + l2^2 - r^2)/ (2*l1*l2));
    q3 = q2 + Alpha - pi/2;
end