
%--------------------------------------------------------------------------
% This script is a) part of mechanical assignment LCSD-03  
% Date - Dec 12 2024 : time : 16:46
% Author - AG Kaimal 
%--------------------------------------------------------------------------
% Declaring all the variables in the question  
%--------------------------------------------------------------------------
% 1.0 Electrical variables  
%syms R real;  % Resistor value Symbolic
R=1; % Resistor value real
Li= 0;  % Inductor value  

% 1.2 Connector variables  
Ke = 0.1;  % Velocity constant  
Kt = 0.1;  % Torque constant  

% 1.3 Mechanical variables  
J1 = 0.00001;  % Moment of inertia of first motor
J2 = 4 * 0.00001;  % Moment of inertia of second motor
%syms D1 real;  % Damping coefficient for first motor Symbolic
D1 = 20;  % Damping coefficient for first motor Symbolic
D2 = 2;   % Damping coefficient for second motor
Bf = 2 * 0.001;  % Friction coefficient

%% Building continous A and b matrices 
%--------------------------------------------------------------------------
%
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

% building C matrix for both the cases 

% 3.1 Case 1: Output y = [theta2; omega2]
C1 = [0 0 1 0 0; 0 0 0 1 0]; 
D1_dir = [0 0;0 0]; % direct input matrix for case-01

% 3.2 Case 2: Output y = [theta2; omega2]
C2 = [0 -Ke/R 0 0 0; 0 0 D2/Bf 0 -D2/Bf];
D2_dir = [-1/R 0; 0 -1/Bf]; % direct input matrix for case-02

%% Descritization of matrix 
Ts=0.001; % sampling time 
A_d=expm(A*Ts);
syms t real
Bfun=@(t)expm(A*t)*B;
B_d=integral(Bfun,0,Ts,'ArrayValued',true);
fprintf('\n\tThis the Disc A matrix: \n');
disp(A_d);
fprintf('\n\tThis the Disc B matrix: \n');
disp(B_d);

%% adding noise to the matrix 
%--------------------------------------------------------------------------
% Noise is added to the state space model 
%--------------------------------------------------------------------------
N=B_d;
fprintf('\n\t This is the N matrix: \n');
disp(N);
%% covarience matrix of noise 
confidence=99.7;
t_bound=0.1;
v_bound=0.3;
den=10*(100-confidence);
volt_var=(v_bound/den)^2;
trq_var=(t_bound/den)^2;
Q=[volt_var 0;0 trq_var];
Cov_NV1=N*Q*pinv(N);          % how the noise maps into state space model
%% discretization of Covarience matrix 
syms tau real
Q_d = integral(@(tau) expm(A * tau) * N * Q * N' * expm(A' * tau), 0, Ts, 'ArrayValued', true);
disp('Discrete-time process noise covariance (Q_d):');
disp(Q_d);

%% measurement noise 
sigma_theta2 = 0.02/3;
sigma_omega2 = 0.01/3;
theta2_noise = sigma_theta2 * randn(1,1);
omega2_noise = sigma_omega2 * randn(1,1);
R= [sigma_theta2^2, 0;
                    0, sigma_omega2^2];
R_d=(1/Ts)*(R);

%% Building Kalman observer 
% basic expression :x_dt=A*x+B*u+N*w
% basic expression :y=C*x+v
% basic expression :x_hat_dt=A*x_hat+B*u+N*w+L*(y-C*x_hat)
%                  :x_til=x-x_hat
%                  :x_til_dt=(A-LC)x_til-Lv
R = eye(size(C1, 1)); % noise covarience matrix for measurement noise 
[P, ~, ~] = dare(A_d', C1',  Q_d , R, [], 'factor'); % Riccati equations 
L = P * C1' * inv(C1 * P * C1' + R_d);                         % Kalman Gain 
fprintf('\n\tThis is the Kalman Gain (L):\n');
disp(L);
observer_eigenvalues = eig(A_d - L * C1);
fprintf('\n\tObserver Eigenvalues:\n');
disp(observer_eigenvalues);
%% LQG controller standard w/o integral actionfor tracking 
% Define cost matrices
Qx_s = eye(size(A_d)); % State weighting
Qu_s = 0.1 * eye(size(B_d, 2)); % Control weighting

% Solve LQR
[P_lqr_s, ~, ~] = dare(A_d, B_d, Qx_s, Qu_s);
K = (B_d' * P_lqr_s * B_d + Qu_s) \ (B_d' * P_lqr_s * A_d);
disp('LQR Gain (K):');
disp(K);
%% LQG controller design with integral action for tracking 
A_a = [A_d, zeros(size(A_d, 1), size(C1, 1));-C1, eye(size(C1, 1))];
B_a = [B_d; zeros(size(C1, 1), size(B_d, 2))];

Qx = blkdiag(eye(size(A_d)), eye(size(C1, 1)));  % State weighting
Qu = 0.1 * eye(size(B_d, 2));                    % Control weighting
% Solve LQR for the augmented system
[P_lqr, ~, ~] = dare(A_a, B_a, Qx, Qu);
K_a = (B_a' * P_lqr * B_a + Qu) \ (B_a' * P_lqr * A_a);
disp('LQR Gain for Augmented System (K_a):');
disp(K_a);