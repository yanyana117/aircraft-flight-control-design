%{
- Xingran Huang
- Aero 625 project
- NZSP
- Sinusoidal Input?
%}

clc;
clear all;
format short ;

H = 100; % feet
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

%% Add actuator dynamics to the overall system

A_new = [A B; zeros(2,5) [-10 0;0 -10]]
% A_new(:,5) = []
B_new = [zeros(5,2); [10 0; 0 10]]
C_new = [diag([1,1,1,1,1,1,1]); zeros(2,7); 0 0 0 0 0 -10 0; 0 0 0 0 0 0 -10]
D_new = [zeros(7,2); 1 0; 0 1; 10 0; 0 10]

% Establish new system 
sys = ss(A_new,B_new,C_new,D_new)

%% Time constant (自己猜)
T = 0.6;

% Set up Q and R
Q = [1e-100        0     0    0       0   0   0; 
         0       1e-100  0    0       0   0   0;
         0         0     1    0       0   0   0;
         0         0     0  1e-100    0   0   0;
         0         0     0    0      30   0   0;
         0         0     0    0       0   3   0;
         0         0     0    0       0   0   5]


R = [1   0; 0   40]

%% Simulate NZSP system
t_final = 20;        % Final time, try 20

table.x = [0; 0; 0; 0; 0; 0; 0]; % Initialize position
table.t = 0: .01: t_final; % Timeline
ym = 5 * (pi/180); % Track 45 degrees


% pacing: Contains the period, T and the value h, the discritization value
pacing.T = T;
pacing.h = 0.01;
h = 0.01;

% Build SDR discrete system and convert to discrete state space model
[sys_SDR,~] = c2d(sys,h)
[gains.K, gains.Q_hat, gains.R_hat, gains.M, gains.S, gains.E] = lqrdjv(A_new,B_new,Q,R,T);

phi = sys_SDR.A;
gamma = sys_SDR.B;
C1 = sys_SDR.C;
D1 = sys_SDR.D;

H = [0 0 0 0 1 0 0];
[pi12, pi22] = QPMCALC(phi-eye(size(phi)), gamma(:,2),H, 0);

% Initialize variables 
frames = 0:T:t_final;

%% Create points for a deflection surface
slope = ym/4; % 1.9 comes from time requirement

table.track = zeros(1,length(table.t)); % Initialize variable

for j = 1:(length(table.t)-1)
    if j < 401
        table.track(1,j+1) = table.track(1,j) + slope*pacing.h;
    else
        table.track(1,j+1) = ym;
    end
end


table.pi12 = pi12
table.pi22 = pi22

% Initial Variables
table.u(:,1) = (pi22 + gains.K*pi12)*ym - gains.K * table.x(:,1); % Update control values
table.y(:,1) = C1 * table.x(:,1) + D1 * table.u(:,1); % Update output values


% Begin Loop over Values:
for i = 1:(length(table.t)-1)
    
    table.x(:,i + 1) = phi * table.x(:,i) + gamma * table.u(:,i); % Update state values
    
    if sum(ismember(table.t(1,i+1),frames)) > 0
        table.u(:,i+1) = (pi22 + gains.K*pi12)*ym - gains.K * table.x(:,i+1); % Update control values
    else
        table.u(:,i+1) = table.u(:,i); % Update control value for ZOH
    end
    
    table.y(:,i+1) = C1 * table.x(:,i+1) + D1 * table.u(:,i+1); % Update output values
end


table.u(:,length(table.t)) = (pi22 + gains.K*pi12)*ym - gains.K * table.x(:,end); % Update control values
table.y(:,length(table.t)) = C1 * table.x(:,end) + D1 * table.u(:,end); % Update output values


%% Characteristics
A_close = A_new - (B_new * gains.K)
[V_close,D_close] = eig(A_close)
sys_NZSP = ss(A_close,B_new,C_new,[0]);
damp(sys_NZSP)

%% plots
Title = 'F-16A Fighting Falcon, lat/d Dynamics';

% 'F-16A Fighting Falcon, lat/d' Dynamics States
figure
hold on
sgtitle(Title + " Dynamics States")

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
% [y_max,index] = max((180/pi)*table.y(10,:)) % highest point in the plot
% xx = table.t(1,:)
% x_max = xx(index)
% text(x_max,y_max,['(',num2str(x_max),',', num2str(y_max),')'],'color','k')     
% plot(table.t(1,:), (180/pi)*table.y(10,:),x_max,y_max,'ro')
plot(table.t(1,:), (180/pi)*table.y(10,:))
ylabel('\delta_a [degs/s]')



subplot(2,1,2)
plot(table.t(1,:),(180/pi)*table.y(11,:))
xlabel('Time [s]')
ylabel('\delta_r [degs/s]')

hold off




% ask do I need to show Rise Time and Settling Time?









