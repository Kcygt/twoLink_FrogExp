% Prefilter (Low-Pass)
Wn = 20;
sigma = 1;

% Desired Joint Trajectory
qDes = [ -0.4240,   2.4189;  
          0.1296    1.9552;  
          0.0,      1.5708;  
         -0.5139    1.9552;  
         -0.4240,   2.4189  
        ]; 


A = [0 0 1 0;
     0 0 0 1;
     -Wn^2 0 -2*sigma*Wn 0;
     0 -Wn^2 0 -2*sigma*Wn ];

u = [qDes(1,1) qDes(1,2) 0 0 ];
