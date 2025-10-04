%--------------------------------------------------------------------------
% This script is the solution for M-2 assignment of 
% Date - Nov 26 2024 
% Author - AG Kaimal 
%--------------------------------------------------------------------------
% Declaring all the variables in the question  
%--------------------------------------------------------------------------
% 1.0 Electrical variables  
syms R real;  % Resistor value Symbolic
L = 0;  % Inductor value  

% 1.2 Connector variables  
Ke = 0.1;  % Velocity constant  
Kt = 0.1;  % Torque constant  

% 1.3 Mechanical variables  
J1 = 0.00001;  % Moment of inertia of first motor
J2 = 4 * 0.00001;  % Moment of inertia of second motor
syms D1 real;  % Damping coefficient for first motor Symbolic
D2 = 2;   % Damping coefficient for second motor
Bf = 2 * 0.001;  % Friction coefficient

%% Building symbolic A and B matrices
%--------------------------------------------------------------------------
% 2.0 Building the A and B matrices  
%--------------------------------------------------------------------------
A = [0 1 0 0 0;
     (-D1/J1) (-Ke*Kt/(J1*R)) (D1/J1) 0 0;
     0 0 0 1 0;
     (D1/J2) 0 -((D2/J2)+(D1/J2)) 0 (D2/J2);
     0 0 (D2/Bf) 0 (-D2/Bf)];

B = [0 0;
     (Kt/(J1*R)) 0;
     0 0;
     0 0;
     0 (1/Bf)];
% Display the A and B matrices
fprintf('\n\tThis is the A matrix:\n'); 
disp(A); 
fprintf('\n\tThis is the B matrix:\n'); 
disp(B);
%% descritization of matrix 

%% Building Controllability matrix 
% controllability matrix 
n=length(A(:,1));   % no. of states in the system

Q = [B A*B (A^2)*B (A^3)*B (A^4)*B];  % general formula 
% Display the Qmatrix 
fprintf('\n\tThis is the Q matrix:\n'); 
disp(Q);

%% - Building Observability Matrix 

% C and D matrix is required for observability 
%--------------------------------------------------------------------------
% 3.0 Building the output C and D matrices 
%--------------------------------------------------------------------------
% 3.1 Case 1: Output y = [theta2; omega2]
C1 = [0 0 1 0 0; 0 0 0 1 0]; 


% 3.2 Case 2: Output y = [theta2; omega2]
C2 = [0 -Ke/R 0 0 0; 0 0 D2/Bf 0 -D2/Bf]; 
D2_dir = [-1/R 0; 0 -1/Bf]; % direct input matrix for case-02  

%--------------------------------------------------------------------------

At=transpose(A);
C1t=transpose(C1);
C2t=transpose(C2);
O1og=[];
O1=[C1t At*C1t At^2*C1t At^3*C1t At^4*C1t]';
O2=[C2t At*C2t At^2*C2t At^3*C2t At^4*C2t]';

% Display the O matrices
fprintf('\n\tThis is the O1 matrices:\n'); 
disp(O1);

fprintf('\n\tThis is the O2 matrices:\n'); 
disp(O2);

%% - native MatLab functions ctrb
An=subs(A,[R D1],[1 20]);
Bn=subs(B,[R D1],[1 20]);
Q_control=ctrb(An,Bn);
fprintf('\n\tThis is the controlability matrix built using native function with values substituted:\n');
disp(Q_control);

C1n=subs(C1,[R D1],[1 20]);
C2n=subs(C2,[R D1],[1 20]);

O1_observe=obsv(An,C1n);
O2_observe=obsv(An,C2n);
fprintf('\n\tThis is the observability matrix built using native function with values substituted:\n');
disp(O1_observe);
fprintf('\n\tThis is the observability matrix built using native function with values substituted:\n');
disp(O2_observe);

%% - Pre loop data
% nested loops for finding values of drop in controllability matrix 

%define the range running R values
Rrun=linspace(0.001,10000,25);

%define the range of running D1 values 
D1run=linspace(0.1,400,25);

% nest loop building 

rankog=rank(double(subs(Q,[R ,D1],[1 ,20])));
ranki=zeros(length(Rrun),length(D1run));
uncntrl_values=[];
safe_value=[];
Ouncntrl_values=[];
Osafe_value=[];
O2uncntrl_values=[];
O2safe_value=[];
%% - Loop 
for i=1:length(Rrun)
    for j=1:length(D1run)
        Qi=subs(Q,[R D1],[Rrun(i) D1run(j)]);
        ranki(i,j)=rank(Qi);
        if ranki(i, j) < rankog
            uncntrl_values=[uncntrl_values; Rrun(i), D1run(j), ranki(i,j)];
        else 
            safe_value=[safe_value; Rrun(i), D1run(j), ranki(i,j)];
        end
    end
end
%% - Loop test for Observability 
for i=1:length(Rrun)
    for j=1:length(D1run)
        Oi=subs(O1,[R D1],[Rrun(i) D1run(j)]);
        ranki(i,j)=rank(Oi);
        if ranki(i, j) < rankog
            Ouncntrl_values=[Ouncntrl_values; Rrun(i), D1run(j), ranki(i,j)];
        else 
            Osafe_value=[Osafe_value; Rrun(i), D1run(j), ranki(i,j)];
        end
    end
end
for i=1:length(Rrun)
    for j=1:length(D1run)
        Oi=subs(O2,[R D1],[Rrun(i) D1run(j)]);
        ranki(i,j)=rank(Oi);
        if ranki(i, j) < rankog
            O2uncntrl_values=[O2uncntrl_values; Rrun(i), D1run(j), ranki(i,j)];
        else 
            O2safe_value=[O2safe_value; Rrun(i), D1run(j), ranki(i,j)];
        end
    end
end

%% Descrete time matrix 
% -------------------------------------------------------------------------
% descritization 
% -------------------------------------------------------------------------

% Building A_d descrete matrix 
% A_d = exp(hA)

% Sampling time interval 
h=0.001;      

A_test=double(subs(A,[R D1],[1 20]));  
B_test=double(subs(B,[R D1],[1 20]));

A_disc_test= expm(h.*A_test);

%--------------------------------------------------------------------------
% Building B_d matrix for zero hold implementation 
%--------------------------------------------------------------------------
syms t real
Bfun=@(t)expm(A_test*t)*B_test;
B_d=integral(Bfun,0,h,'ArrayValued',true);

%--------------------------------------------------------------------------
% Stability analysis 
%--------------------------------------------------------------------------

A_d1=double(subs(A,[R,D1],[1,20]));
Stab_eig=eig(A_disc_test);
fprintf('\n\tThe eigen values of discretized A matrix with substituting values: \n');
disp(Stab_eig);
%%
%--------------------------------------------------------------------------
% Contrallability matrix of the discrete system
%--------------------------------------------------------------------------
A_dt=A_disc_test;
Q_d=[B_d A_dt*B_d A_dt^2*B_d A_dt^3*B_d A_dt^4*B_d];

%--------------------------------------------------------------------------
