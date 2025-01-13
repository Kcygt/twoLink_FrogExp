%% Two link with prefilter
% Input is Uinp, prefilter states are [u1p u2p u1 u2 q1p q2p q1 q2]'.
%
% use help twolink with prefilter to get demo commands
 % x0=zeros(8,1);
 % x0(1:2)=[.2;-.1];
 % [t,y]=ode113(@(t,x) twolinkwithprefilter(t,x,1),[0 7],x0);
 % figure(1);plot(t,y(:,1:4));shg;  figure(2);plot(t,y(:,5:8));shg

 
% >> [t,y]=ode45(@(t,x) twolinkwithprefilter(t,x,1),[0 17],x0);


function dxdt=twolinkwithprefilter(t,x,wn)
% Prefilter input is assumed to be added to Uinp(3:4) (should probably
% do it via a B vector so it is more similar to xdot=Ax+Bu)
%
% Here robot link lengths are in this function but could be passed in
% As of 13/1/2025 it is untested and  unstable

zeta=1;
A=[zeros([2 2]) eye(2);-eye(2)*wn^2 -eye(2)*2*zeta*wn]; % note the wn^2 !!

%% Set up a two link arm
q=x(5:6);
qd=x(7:8);
q1p=x(7);q2p=x(8);
q1=x(5);q2=x(6);

bj=.01; % energy dissipation term
kj=.3; % controller gain
L_1=1;
L_2=1;
m_1=.5;
m_2=.5;
% derived constants
ka=L_2^2*m_2;
kb=(L_1*L_2*m_2);
kc=L_1^2*(m_1+m_2);

M=[ka+2*kb*cos(q2)+kc  ka+kb*cos(q2);
   ka+kb*cos(q2) ka];
V=ka*sin(q2)*([0 -1;1 0]*[q1p^2;q2p^2]+[-2*q1p*q2p;0]);

Numerator=V+[bj 0;0 bj]*qd + [kj 0;0 kj]*(q-x(1:2));
%Numerator=V;
qdd=M\Numerator;

if t<.4;
    dotx=[0;0;0;0];
elseif (t<5)
    Uinp=[0;0;1;1.5]*wn^2; 
    dotx=A*x(1:4)+Uinp;
else
    Uinp=[0;0;0.2;0]*wn^2;
    dotx=A*x(1:4)+Uinp;
end


dxdt=[dotx;qd;qdd];


end

