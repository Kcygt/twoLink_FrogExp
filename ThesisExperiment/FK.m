function [x,y,z] = FK(q1,q2,q3)
    l1 = 0.208; 
    l2 = 0.168;  
    x = sin(q1).*(l1*cos(q2)+l2*sin(q3));
    y = l2-l2*cos(q3)+l1*sin(q2);
    z = -l1+cos(q1).*(l1*cos(q2)+l2*sin(q3));
end

