 [fX,fY,cartesianX,cartesianY] = defineSquare(0.9, 0.9, 1);
 qDes = inverse_kinematics(cartesianX,cartesianY, 1, 1)';
 time = [5, 10, 15, 20, 30];
 x0=zeros(8,1);
 x0(1:2)=[qDes(1,1);qDes(1,2)];
 
 wn = 2;
 
 [t,y]=ode113(@(t,x) myTwolinkwithprefilter(t,x,wn,time, qDes),[0 30],x0);
 % figure(1);plot(t,y(:,1:4));shg;  figure(2);plot(t,y(:,5:8));shg
 xAct =  forward_kinematics(y(:, 5), y(:, 6), 1, 1);
 xDes =  forward_kinematics(qDes(:, 1), qDes(:, 2), 1, 1);
 figure(3); plot(xAct(:,1),xAct(:,2),'-'); hold on; plot(xDes(:,1),xDes(:,2),'o-')