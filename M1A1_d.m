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

% Display the A and B matrices
fprintf('\n\tThis is the A matrix:\n'); 
disp(A); 
fprintf('\n\tThis is the B matrix:\n'); 
disp(B);

% Calculate the eigenvalues of matrix A
A_eigen = eig(A); 
fprintf('\nThese are the eigenvalues of A:\n'); 
disp(A_eigen);

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


%--------------------------------------------------------------------------
% input output stability | poles zeros cancelation 
%--------------------------------------------------------------------------
p1=pole(sys1);
z1=tzero(sys1(1,1));z2=tzero(sys1(1,2));z3=tzero(sys1(2,1));z4=tzero(sys1(2,2));

%combination of all the zeros 

z1e=[z1;z2;z3;z4];z1e=unique(z1e,"rows");

p2=pole(sys2);
z1_2=tzero(sys2(1,1));z2_2=tzero(sys2(1,2));z3_2=tzero(sys2(2,1));z4_2=tzero(sys2(2,2));

%combination of all the zeros 

z2e=[z1_2;z2_2;z3_2;z4_2];z2e=unique(z2e,"rows");

%--------------------------------------------------------------------------
% pole-zero calcellation for case 1
%--------------------------------------------------------------------------
tol = 0.000001;
count1=0;
for i = 1:length(z1e)
    for j = 1:length(p1)
        if abs(p1(j)-z1e(i))<tol
            count1=count1+1;
        end
    end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% pole-zero calcellation for case 2
%--------------------------------------------------------------------------
tol2 = 0.000001;
count2=0;
for i = 1:length(z2e)
    for j = 1:length(p2)
        if abs(p2(j)-z2e(i))<tol2
            count2=count2+1;
        end
    end
end
%--------------------------------------------------------------------------
