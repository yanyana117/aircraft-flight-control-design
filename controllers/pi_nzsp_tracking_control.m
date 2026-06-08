clc
clear all

clear table;
%{
- Xingran Huang
- Aero 625 project
- PI-NZSP
%}

%% Given initail continuous system (A,B)
A = [ -0.132    0.324    -0.94      0.149  0;...
      -10.614  -1.179     1.0023    0      0; ...
       0.997   -0.00182  -0.259     0      0;...
       0         1        0.34      0      0;...
       0         0        1.0561    0      0]; 
% A B for actuators：make somthing move or operate
B = [0.0069 0.0189; -5.935  1.203; -0.122   -0.614; 0    0; 0    0];
C = diag([1,1,1,1,1]);
D = zeros(5,2); % 5 rwo, 2 col all zero
% H = 100; % feet
H = [0 0 0 0 1 0 0]; % Define the H Matrix

%% Add actuator dynamics to the overall system
A_new = [A B; zeros(2,5) [-10 0;0 -10]];
B_new = [zeros(5,2); [10 0; 0 10]];
C_new = [diag([1,1,1,1,1,1,1]); zeros(2,7); 0 0 0 0 0 -10 0; 0 0 0 0 0 0 -10];
D_new = [zeros(7,2); 1 0; 0 1; 10 0; 0 10];

% Establish new system 
sys = ss(A_new,B_new,C_new,D_new);

%% Q & R & Time constant selected by self
% smaller intrgal gains:Q
Q = diag([1e-100, 0.001, 200, 0.1, 0.1, 0.1, 0.1, 1]); 
R = [9 0; 0 70];
T = 0.01;
t_final = 10;        % Final time, try 20
table.x = [0; 0; 0; 0; 0; 0; 0; 0]; % Initialize position
table.t = 0: .01: t_final; % Timeline
ym = 30 * (pi/180); % cheange to Track 30 degrees
Disturbance = [0;0.01;0;0;0;0;0;0];
pacing.T = 0.05; % smaller
% pacing: Contains the period, T and the value h, the discritization value
% pacing.T = T;
pacing.h = 0.01;
h = 0.01;

%% Simulate PI-NZSP system
N = length(sys.A);

%% Build discrete system for PI matrix: 
% Solves for the gains and discrete system of the continuous system

% Set up discrete system for PI
[n,~] = size(sys.A);
[~,b] = size(sys.B);
[c,~] = size(H);

% Continuous System
% A = sys.A;    % phi
% B = sys.B;  % gamma
% C = sys.C;     % C1
% D = sys.D;     % D1

A1 = [sys.A, zeros(n,c); H, zeros(c,c)]
B1 = [sys.B; zeros(c,b)]

% Discrete
[sys_PI_NZSP,~] = c2d(sys,h);

%% PI Matrix
% PI added to discrete system
% sys_PI_NZSP_discrete = 

[n,~] = size(sys_PI_NZSP.A);
[~,b] = size(sys_PI_NZSP.B);
[c,~] = size(H);
[e,f] = size(sys_PI_NZSP.C);
[~,k] = size(sys_PI_NZSP.D);

% 这组pi是pi matrix里的ABCD
A_pi = [sys_PI_NZSP.A zeros(n,c); (H * T) eye(c)];
B_pi = [sys_PI_NZSP.B; zeros(c,b)];
C_pi = [sys_PI_NZSP.C zeros(e,c); zeros(c,(f+c))];
C_pi(e+c,f+c) = 1;
D_pi = [sys_PI_NZSP.D; zeros(c,k)]; 

% Contains the discrete version of the continous system
% PI added to discrete system
% sys_PI是 pi_matrix里的sys2 (output)，
%          build_dis_sys_pi里的sys3 (output)

sys_combine = ss(A_pi,B_pi,C_pi,D_pi);
% gaiming pi matrix output

[gains.K, gains.Q_hat, gains.R_hat, gains.M, gains.S, gains.E] = lqrdjv(A1,B1,Q,R,pacing.T);


% New sys for PI-NZSP (pi_nzsp 文件)
A_new_sys = sys_combine.A;
B_new_sys = sys_combine.B;
C_new_sys = sys_combine.C;
D_new_sys = sys_combine.D;
H_new_sys = [0 0 0 0 1 0 0];

% Old way for doing NZSP
A_old = sys_PI_NZSP.A;
B_old = sys_PI_NZSP.B;
C_old = sys_PI_NZSP.C;
D_old = sys_PI_NZSP.D;

T = pacing.T;
h = pacing.h;

Y = [zeros(length(A_new_sys)-1,1);-h]; % Goes into x equation

% Initialize variables 
frames = 0:T:t_final;

% Create points for a deflection surface
slope = ym/4; % 1.9 comes from time requirement

table.track = zeros(1,length(table.t)); % Initialize variable
for j = 1:(length(table.t)-1)
    if j < 401
        table.track(1,j+1) = table.track(1,j) + slope*pacing.h;
    else
        table.track(1,j+1) = ym;
    end
end



%% Solve for pi22 and pi12 
% Check: A_old is phi and B_old is gamma

% [sys_PI_NZSP,gamma_B] = c2d(sys,h);
[pi12,pi22] = QPMCALC(A_old - eye(size(A_old)),B_old(:,2),H,0);

% [pi12,pi22] = QPMCALC(sys_PI_NZSP - eye(size(sys_PI_NZSP)),gamma_B(:,2),H,0);

table.pi12 = pi12;
table.pi22 = pi22;


% Initial Variables
kkk = gains.K(:,(1:N));
table.u(:,1) = (pi22 + kkk * pi12)*ym - gains.K * table.x(:,1);  % Update control values
table.y(:,1) = C_new_sys * table.x(:,1) + D_new_sys * table.u(:,1); % Update output values


% Begin Loop over Values
for i = 1:(length(table.t)-1)
    
    table.x(:,i + 1) = A_new_sys * table.x(:,i) + B_new_sys * table.u(:,i) + ...
        (Y*ym) + Disturbance * h; % Update state values
    
    if sum(ismember(table.t(1,i+1),frames)) > 0
        table.u(:,i+1) = (pi22 + gains.K(:,(1:N))*pi12)*ym - gains.K * table.x(:,i+1); % Update control values
    else
        table.u(:,i+1) = table.u(:,i); % Update control value for ZOH
    end
    
    table.y(:,i+1) = C_new_sys * table.x(:,i+1) + D_new_sys * table.u(:,i+1); % Update output values
end


% Update the last values
table.u(:,length(table.t)) = (pi22 + gains.K(:,(1:N))*pi12)*ym - gains.K * table.x(:,end); % Update control values
table.y(:,length(table.t)) = C_new_sys * table.x(:,end) + D_new_sys * table.u(:,end); % Update output values



%% plots
Title = 'F-16A Fighting Falcon, lat/d Dynamics';

% 'F-16A Fighting Falcon, lat/d' Dynamics States
figure
hold on
sgtitle(Title + "  States")

subplot(5,1,1)
plot(table.t(1,:),(180/pi)*table.y(1,:))
ylabel('State \beta [degs]')

subplot(5,1,2)
plot(table.t(1,:),(180/pi)*table.y(2,:))
ylabel('State p [degs/s]')

subplot(5,1,3)
plot(table.t(1,:),(180/pi)*table.y(3,:))
ylabel('State r [degs/s]')

subplot(5,1,4)
plot(table.t(1,:),(180/pi)*table.y(4,:))
ylabel('State \phi [degs]')

subplot(5,1,5)
hold on
plot(table.t(1,:),(180/pi)*table.y(5,:))
plot(table.t(1,:),(180/pi)*table.track(1,:))
legend('State','Desired Step Output')
hold off

xlabel('Time [s]')
ylabel('State \psi [degs]')
hold off

%% 'F-16A Fighting Falcon, lat/d' Control Graphs
% Plot the Control Graphs
figure
hold on
sgtitle(Title + " Controls")

subplot(2,1,1)
plot(table.t(1,:),(180/pi)*table.y(6,:))
ylabel('\delta_a [degs]')

subplot(2,1,2)
plot(table.t(1,:),(180/pi)*table.y(7,:))

xlabel('Time [s]')
ylabel('\delta_r [degs]')
hold off

%% 'F-16A Fighting Falcon, lat/d' Commands Graphs
figure
hold on
sgtitle(Title + " Commands")

subplot(2,1,1)
plot(table.t(1,:),(180/pi)*table.y(8,:))
ylabel('\delta_a [degs]')

subplot(2,1,2)
plot(table.t(1,:),(180/pi)*table.y(9,:))

xlabel('Time [s]')
ylabel('\delta_r [degs]')
hold off

%% 'F-16A Fighting Falcon, lat/d'  Rates Graphs
figure
hold on
sgtitle(Title + " Rates")

subplot(2,1,1)
plot(table.t(1,:), (180/pi)*table.y(10,:))
ylabel('\delta_a [degs/s]')

subplot(2,1,2)
plot(table.t(1,:),(180/pi)*table.y(11,:))
xlabel('Time [s]')
ylabel('\delta_r [degs/s]')

hold off

%%

Q_hat = gains.Q_hat;
R_hat1 = gains.R_hat;
M1 = gains.M;
K = gains.K

% Steady State Values
x_star = table.pi12 * ym;
u_star = table.pi22 * ym;

% Closed Loop Matrix Characteristics
A_cl = sys_combine.A - (sys_combine.B * K);
[V_cl,D_cl] = eig(A_cl)
sys_final = ss(A_cl,sys_combine.B,sys_combine.C,sys_combine.D)
damp(sys_final)





