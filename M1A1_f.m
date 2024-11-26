
%--------------------------------------------------------------------------
% This script is d), e), f) part of mechanical assignment LCSD  
% Date - Nov 18 2024 
% Author - AG Kaimal 
%--------------------------------------------------------------------------
% Declaring all the variables in the question  
%--------------------------------------------------------------------------
% 1.0 Electrical variables  
R = 1;  % Resistor value
L = 0;  % Inductor value  

% 1.2 Connector variables  
Ke = 0.1;  % Velocity constant  
Kt = 0.1;  % Torque constant  

% 1.3 Mechanical variables  
J1 = 0.00001;  % Moment of inertia of first motor
J2 = 4 * 0.00001;  % Moment of inertia of second motor
D1 = 20;  % Damping coefficient for first motor
D2 = 2;   % Damping coefficient for second motor
Bf = 2 * 0.001;  % Friction coefficient  

%--------------------------------------------------------------------------
% 2.0 Building the A and B matrices  
%--------------------------------------------------------------------------
A = [0 1 0 0 0;
     (-D1/J1) (-Ke*Kt/(J1*R)) (D1/J1) 0 0;
     0 0 0 1 0;
     (D1/J2) 0 -((D2/J2)+(D1/J2)) 0 (D2/J2);
     0 0 (D2/Bf) 0 (-D2/Bf)];

B = [0 0;
     (-Kt/(J1*R)) 0;
     0 0;
     0 0;
     0 (-1/Bf)];

%--------------------------------------------------------------------------
% 3.0 Building the output C and D matrices 
%--------------------------------------------------------------------------
% 3.1 Case 1: Output y = [theta2; omega2]
C1d = [0 0 1 0 0; 0 0 0 1 0]; 
D1d = [0 0; 0 0]; 

% 3.2 Case 2: Output y = [theta2; omega2]
C2d = [0 -Ke/R 0 0 0; 0 0 D2/Bf 0 -D2/Bf]; 
D2d = [-1/R 0; 0 -1/Bf]; 

%--------------------------------------------------------------------------
% building State Space model 
%--------------------------------------------------------------------------
sys1=ss(A,B,C1d,D1d);
sys2=ss(A,B,C2d,D2d);
    
t = 0:0.01:2;
u = zeros(size(t),2);
x0 = [0.01 0 0];
[y,t,x] = lsim(sys1,u,t,x0);
plot(t,y)