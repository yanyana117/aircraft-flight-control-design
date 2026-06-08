%{
- Xingran Huang
- Aero 625 project

- PI_SDR controller
- Sinusoidal Input?

By consider the lateral/directional F-16A Fight Falcon Linear Model
%}

clc
clear all
format short


%% Initial Condition:
M1 = 0.18;
U1 = 205;   % feet/sec
H1 = 100;   % feet
alhpha_1 = 18.8;    % degree
q_bar = 50; % psf
% x_cg = 0.35*c_bar;
T = 0.05;   % sample time

% Given initail continuous system (A,B)
A = [ -0.132    0.324    -0.94      0.149  0;...
      -10.614  -1.179     1.0023    0      0; ...
       0.997   -0.00182  -0.259     0      0;...
       0         1        0.34      0      0;...
       0         0        1.0561    0      0]; 
% A B for actuators：make somthing move or operate
B = [0.0069 0.0189; -5.935  1.203; -0.122   -0.614; 0    0; 0    0];

%% Open Loop Characteristics
% Eigenvalues and eigenvectors
[V_open,D_open] = eig(A);

% lambda — Eigenvalues (returned as vector) column vector
lambda = eig(A);
damp(A)

% Time constant
tau = -1 ./ nonzeros(real(lambda));

% wn:natural freq; zeta: damping ratio; p: poles
[wn, zeta, p] = damp(A);

% Sampling frequency: ws > 2*max(wn)
ws = 10*max(wn);
h = .01;

% Sample time
ts = 2 * pi / ws;

%% Add actuator dynamics
% A & B Matrix with actuator dynamics
A_new = [A B; zeros(2,5) [-10 0;0 -10]];
B_new = [zeros(5,2); [10 0; 0 10]];
C_new = [diag([1,1,1,1,1,1,1]); zeros(2,7); 0 0 0 0 0 -10 0; 0 0 0 0 0 0 -10];
D_new = [zeros(7,2); 1 0; 0 1; 10 0; 0 10];

% Establish new system 
sys = ss(A_new,B_new,C_new,D_new);
% 7b: Controller gains and the weighting matrices used to calculate them.

% How to get Q & R? Define by us,tsting start from diag([1,1,1,1,1,1])
% Designer selects continuous Q and R, algorithm calculates dQ，dR, M.

Q = diag([0.4, 0.1,1, 10, 100, 1, 1, 350])
R = [1250 0; 0 70]

CC0 = [0;0;0;0;0;0;0];
C1 = [diag([1, 1, 1, 1, 1, 1, 1]) CC0; 0 0 0 0 0 -10 0 0; ...
      0 0 0 0 0 0 -10 0; zeros(2,8); 0 0 0 0 0 0 0 1];  

D1 = [zeros(7,2); 10 0; 0 10; 1 0; 0 1; 0 0];

%% Add integrated state to A, and B, C, D
A = A_new;
B = B_new;
C = C_new;
D = D_new;
A_y = [0 0 0 0 1 0 0 0 ];
A_int = [A, [0 0 0 0 0 0 0 ]' ; A_y];

B_y = [0 0];
B_int = [B; B_y];

%% Setup discretized simulation
y_m = 0; %desired position (regulator problem)
x_k = [0 ; 0 ; 0 ; -12 ; 0; 0; 0; 0]; % initial conditions
x_k(2:7) = x_k(2:7)*(pi/180); % convert to radians for simulation

T = .01;
to = 0;  % initial time
tf = 10; % final time
h = .01;

% Discretized simulation
[phi, gamma] = c2d(A, B, h);

%% Add integrated states to Phi Gamma
phi_y = [0 0 0 0 T 0 0 1];
phi_int = [phi, [0 0 0 0 0 0 0]'; phi_y];
gamma_y = [0 0];
gamma_int = [gamma; gamma_y]

%% Check Controllability / Reachability
% if r = n then A matrix is fully controllable.
% If A is fully controllable then it is also fully reachable

% Discretized controllability matrix
CO = ctrb(phi_int, gamma_int); 
r = rank(CO); 
n = size(A); 

%% Check Observability
% rank = n therefore the system is fully observable and fully constructable
C = eye(8); 
O = obsv(phi_int,C);
r_o = rank(O); 

%%  SDR Cost Function Minimization to Find Gain values for K
nn = eye(8,2)*0;
[k,Qd,Rd,Nd,s,e] = lqrdjv(A_int,B_int,Q,R,nn,T); 
u_k = y_m - (k *x_k);

%% Initialize Variables at t_o
hold_count = 0;
countmax = T/h;
frames = (tf - to)/h;
t = to;
y = (C1 * x_k + D1 * u_k)';

%% Modal Characteristics for closed loop continuous eigen 
damp(A_int - B_int*k) 
[V_close, D_close] = eig(A_int - B_int*k)

%% disturbance
w = .01;
d = [0 0 0 0 T 0 0 0 ]';

%% Run Simulation
for i=1:frames
    if hold_count == countmax
        u_k = y_m - (k*x_k);
        hold_count = 0;
    end

    x_k = phi_int * x_k + gamma_int * u_k +d*w;
    y = [y; (C1*x_k + D1*u_k)'];
    t = [t; i*h];
    hold_count = hold_count + 1;

end


%% Convert output back to degrees
y(:,2:11) = y(:,2:11)*(180/pi);

%% plot
% lateral/directional F-16A Fight Falcon Linear Model

Title = 'F-16A Fighting Falcon, lat/d Dynamics';

figure(1)
hold on
sgtitle(Title + "  States")

subplot(5,2,1)
plot(t,y(:,1))
title('Sideslip Rate Dynamics')
xlabel('time')
ylabel('Sideslip angle \beta [degs]')


subplot(5,2,3)
plot(t,y(:,2))
xlabel('time')
title('Roll Rate Dynamics')
ylabel('Roll Rate p [deg/s]')

subplot(5,2,5)
plot(t,y(:,3))
title('Yaw Rate Dynamics')
xlabel('time')
ylabel('Yaw Rate r [deg/s]')

subplot(5,2,7)
plot(t,y(:,4))
title('Roll Angle Dynamics')
xlabel('time')
ylabel('Roll Angle ϕ [deg]')


subplot(5,2,9)
plot(t,y(:,5))
title('Yaw Angle Dynamics')
xlabel('time')
ylabel('Yaw Angle Ψ [deg]')
hold off


subplot(5,2,2);
plot(t,y(:,6))
title('Actuator Control')
xlabel('time')
ylabel('Actuator Angle [deg]')
grid on;
hold on;
plot(t,y(:,7))
xlabel('time')
legend("Aileron Control" , "Rudder Control")


subplot(5,2,4);
plot(t,y(:,8))
title('Actuator Control Rate')
xlabel('time')
ylabel('Actuator Rate [deg/s]')
grid on;
hold on;
plot(t,y(:,9))
xlabel('time')
legend("Aileron Control Rate" , "Rudder Control Rate")


subplot(5,2,6);
plot(t,y(:,10))
title('Actuator Commanded Value')
xlabel('time')
ylabel('Actuator Command [deg]')
grid on;
hold on;
plot(t,y(:,11))
xlabel('time')
legend("Commanded Aileron Angle" , "Commanded Rudder Angle")


subplot(5,2,8)
plot(t,y(:,1))
title('Integrated State')
ylabel('Y_I')
xlabel('time')
grid on;
hold on;


grid on;
hold on;


