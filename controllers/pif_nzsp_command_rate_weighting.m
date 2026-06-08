%{
- Xingran Huang
- Aero 625 project

- PIF_NZSP_CRW controller
- Sinusoidal Input?

By consider the lateral/directional F-16A Fight Falcon Linear Model
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
A_act = [A B; zeros(2,5) [-10 0;0 -10]]
B_act = [zeros(5,2); [10 0; 0 10]]
C_act = [diag([1,1,1,1,1,1,1]); zeros(2,7); 0 0 0 0 0 -10 0; 0 0 0 0 0 0 -10]
D_act = [zeros(7,2); 1 0; 0 1; 10 0; 0 10]

% Establish new system 
sys = ss(A_act,B_act,C_act,D_act)

%% Q & R & Time constant selected by self
% smaller intrgal gains:Q
% Q = diag([1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1]); 
Q =     [1  0  0   0  0    0  0   0  0  0; 
         0  1  0   0  0    0  0   0  0  0;
         0  0  15  0  0    0  0   0  0  0;
         0  0  0   1  0    0  0   0  0  0;
         0  0  0   0  100  0  0   0  0  0; %end of states
         0  0  0   0  0   10  0   0  0  0; 
         0  0  0   0  0    0  10  0  0  0;
         0  0  0   0  0    0  0  20  0  0;
         0  0  0   0  0    0  0   0  1  0;
         0  0  0   0  0    0  0   0  0  1]

R = [9 0; 0 70]
% R = [1  0; 0 3.2];

T = 0.1;
t_final = 20;        % Final time, try 20
% table.x = [0; 0; 0; 0; 0; 0; 0; 0]; % Initialize position
table.x = [0;0;0;0;0;0;0;0;0;0]; % Initialize position
table.t = 0: .01: t_final; % Timeline
ym = 5 * (pi/180); % Track 5 degrees

Disturbance = [0;0.5;0;0;0;0;0;0;0;0];
pacing.T = 0.05; % smaller
% pacing: Contains the period, T and the value h, the discritization value
% pacing.T = T;
pacing.h = 0.01;
h = 0.01;



    %% Simulate System
    % PIF - NZSP - CRW system Set-Up
    
    N = length(A_act);

        %% build 文件
        [n,~] = size(A_act);
        [~,b] = size(B_act);
        [c,~] = size(H);
        
        
        % Continuous System
        A_lqr = [A_act, B_act, zeros(n,c);  zeros(b,n+b+c); H, zeros(c,b+c)];
        B_lqr = [zeros(n,b);eye(b); zeros(c,b)];
        
        
        
        % Discrete
        [sys_old,~] = c2d(sys,h);
        
                %% PI Matrix
                        
                % Build Systems
                [n,~] = size(sys_old.A);
                [~,b] = size(sys_old.B);
                [c,~] = size(H);
                [e,f] = size(sys_old.C);
                [~,k] = size(sys_old.D);
                
                A_PI_old = [sys_old.A,     sys_old.B,   zeros(n,c); 
                zeros(b,n), eye(b),  zeros(b,c);
                T*H,     zeros(c,b), eye(c)];
                
                B_PI_old = [zeros(n,b);
                T*eye(b);
                zeros(c,b)];
                
                
                C_PI_old = [eye(n+b+c)];
                D_PI_old = [zeros(c+b+n,k)]; 
                sys_PI = ss(A_PI_old,B_PI_old,C_PI_old,D_PI_old);
                sys_PI.Ts = T;

        %% 继续build
        [sys_PI_d2c] = d2c(sys_PI);
        [gains.K, gains.Q_hat, gains.R_hat, gains.M, gains.S, gains.E] = lqrdjv(A_lqr,B_lqr,Q,R,T);


    %% 继续 PIF
    % New System for PI NZSP
    A_PI = sys_PI.A;
    B_PI = sys_PI.B;
    C_PI = sys_PI.C;
    D_PI = sys_PI.D;

    % H = sys.H;
    
    % From old way of doing NZSP
    A1 = sys_old.A;
    B1 = sys_old.B;
    C1 = sys_old.C;
    D1 = sys_old.D;
    
    Y = [zeros(length(A_PI)-1,1);h]; % Goes into x equation
    
    % Initialize variables 
    frames = 0:T:t_final;
    
    % Create points for a deflection surface
    slope = ym/2.5; % 1.9 comes from time requirement
    
    table.track = zeros(1,length(table.t)); % Initialize variable
    
    for j = 1:(length(table.t)-1)
        if j < 250
            table.track(1,j+1) = table.track(1,j) + slope*pacing.h;
        else
            table.track(1,j+1) = ym;
        end
    end
        
    
    % Solve for pi22 and pi12
    [pi12, pi22] = QPMCALC(A1 - eye(size(A1)),B1(:,2),H,0);
    table.pi12 = pi12;
    table.pi22 = pi22;
    
    
    % Initial Variables
    table.u(:,1) = (pi22*gains.K(:,N+1) + gains.K(:,(1:N))*pi12)*ym - gains.K * table.x(:,1); % Update control values
    table.y(:,1) = C_PI * table.x(:,1) + D_PI * table.u(:,1); % Update output values
    
    
    % Begin Loop over Values
    for i = 1:(length(table.t)-1)
        
        table.x(:,i + 1) = A_PI * table.x(:,i) + B_PI * table.u(:,i) - (Y*ym) + Disturbance * h; % Update state values
        
        if sum(ismember(table.t(1,i+1),frames)) > 0
            table.u(:,i+1) = (pi22*gains.K(:,N+1) + gains.K(:,(1:N))*pi12)*ym - gains.K * table.x(:,i+1); % Update control values
        else
            table.u(:,i+1) = table.u(:,i); % Update control value for ZOH
        end
        
        table.y(:,i+1) = C_PI * table.x(:,i+1) + D_PI * table.u(:,i+1); % Update output values
    end
    
    % Update the last values
    table.u(:,length(table.t)) = (pi22*gains.K(:,N+1) + gains.K(:,(1:N))*pi12)*ym - gains.K * table.x(:,end); % Update control values
    table.y(:,length(table.t)) = C_PI * table.x(:,end) + D_PI * table.u(:,end); % Update output values

    %% Simulate System   
Q_hat = gains.Q_hat
R_hat1 = gains.R_hat
M1 = gains.M;
K = gains.K

% Steady State Values
x_star = table.pi12 * ym
u_star = table.pi22 * ym

% New State Equations
A_new = sys_PI.A;
B_new = sys_PI.B;

% Solve for Closed Loop Matrix
A_cl = sys_PI.A - (sys_PI.B * K);

% Eigenvalues (D) & Eigenvectors (V)
[V_cl,D_cl] = eig(A_cl)

% Natural Frequencies & Damping Ratios
sys_final = ss(A_cl,sys_PI.B,sys_PI.C,sys_PI.D)
damp(sys_final) 

%% Plot
  
Title = "F-16A Fight Falcon Lateral Dynamics"';


    figure
    hold on
    sgtitle(Title + " States")
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
    
    % Plot the Commands Graphs
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
    
    % Plot the Rates Graphs
    figure
    hold on
    sgtitle(Title + " Rates")
    subplot(2,1,1)
    plot(table.t(1,:),(180/pi)*table.u(1,:))
    ylabel('\delta_a [degs/s]')
    subplot(2,1,2)
    plot(table.t(1,:),(180/pi)*table.u(2,:))
    xlabel('Time [s]')
    ylabel('\delta_r [degs/s]')
    hold off









