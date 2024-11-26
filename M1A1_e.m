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
% Part-e: Transfer function for Ae and Be
%--------------------------------------------------------------------------
Ae = [0 1 0 0 0;
      (-D1/J1) (-Ke*Kt/(J1*R)) (D1/J1) 0 0;
      0 0 0 1 0;
      (D1/J2) 0 -((D1/J2) + (D2/J2)) 0 (D2/J2);
      0 0 (D2/Bf) 0 (-D2/Bf)];

Be = [0 0;
      (-Kt/(J1*R)) 0;
      0 0;
      0 0;
      0 (-0/Bf)];

%--------------------------------------------------------------------------
% Display the Ae and Be matrices
%--------------------------------------------------------------------------
fprintf('\n\tThis is the Ae matrix:\n'); 
disp(Ae); 

fprintf('\n\tThis is the Be matrix:\n'); 
disp(Be); 

%--------------------------------------------------------------------------
% Calculate the eigenvalues of Ae
%--------------------------------------------------------------------------
Ae_eigen = eig(Ae); 
fprintf('\nThese are the eigenvalues of Ae:\n'); 
disp(Ae_eigen); 

%--------------------------------------------------------------------------
%  Building the output C and D matrices 
%--------------------------------------------------------------------------
% 3.1 Case 1: Output y = [theta2; omega2]
C1e = [0 0 1 0 0; 0 0 0 1 0]; 
D1e = [0 0; 0 0]; 

% 3.2 Case 2: Output y = [theta2; omega2]
C2e = [0 -Ke/R 0 0 0; 0 0 D2/Bf 0 -D2/Bf]; 
D2e = [-1/R 0; 0 0]; 


%--------------------------------------------------------------------------
syse2=ss(Ae,Be,C2e,D2e);
tf(syse2)

%--------------------------------------------------------------------------
% input output stability | poles zeros cancelation 
%--------------------------------------------------------------------------
p1e=pole(syse2);
z1e=tzero(syse2(1,1));z2e=tzero(syse2(1,2));z3e=tzero(syse2(2,1));z4e=tzero(syse2(2,2));

%combination of all the zeros 

z1e_supr=[z1e;z2e;z3e;z4e];z1e_supr=unique(z1e_supr,"rows");

%--------------------------------------------------------------------------