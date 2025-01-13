%% Demo of a critically damped second order filter
%
% Input is Uinp, prefilter states are [u1p u2p u1 u2]'. An easy
% variation of states would be [u1p u1 u2p u2]'. An easy
%
% wn=2;
% [t,x]=ode45(@(t,x) prefilterplay(t,x,wn),[0 30],[0;0;0;0]);plot(t,x);shg;grid on
%
% Probably easiest to set wn high to get close to no prefilter
%
%[t,x]=ode45(@(t,x) prefilterplay(t,x,wn),[0 30],[0;0;0;0]);plot(t,x(:,1:2));shg;grid on
% hold on
% wn=4;
% [t,x]=ode45(@(t,x) prefilterplay(t,x,wn),[0 30],[0;0;0;0]);plot(t,x);shg;grid on
% wn=2.5;
% [t,x]=ode45(@(t,x) prefilterplay(t,x,wn),[0 30],[0;0;0;0]);plot(t,x);shg;grid on
% hold off
%
% the state function needs to encode the change of inputs with time.
% It is probable that some ode solvers might allow an input function to be called
% to do this possibly via the opts function. Here we do it with
% 'if-elseif - ...' statements
function dotx=prefilterdemo(t,x,wn)
% assume the function is f/m= \ddot{x}+b/m \dot{x} = k/m(x-x0)
% then if f=0 we get \ddot{x}=wn^2 x0 -wn^2 x -2\zeta wn \dot{x}

zeta=1;
A=[zeros([2 2]) eye(2);-eye(2)*wn^2 -eye(2)*2*zeta*wn]; % note the wn^2 !!

if t<.5;
    dotx=[0;0;0;0];
elseif (t<5)
    Uinp=[0;0;1;1.5]*wn^2; 
    dotx=A*x+Uinp;
else
    Uinp=[0;0;0.2;0]*wn^2;
    dotx=A*x+Uinp;
end

end

