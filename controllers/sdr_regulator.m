%{
- Xingran Huang
- Aero 625 project

- Introduction From hw 4 problem 7:
For the Project, implement parts 2a, 3, and the 4.1 specifications for the regulator problem only. Please
refer to the project handout. The minimum items you are to provide are:
(a) Closed-loop characteristics:
    eigenvalues
    eigenvectors
    natural frequencies, damping ratios and time constants as applicable
(b) Controller gains and the weighting matrices used to calculate them
(c) Simulation time histories with the digital controller active, showing that you satisfied the specifications,
and including:
    Plots of all commands
    Plots of all state variables
    Plots of all control variables

2. Methodology: For the same plant, design a discrete controller of each type listed below:
a) SDR

3. Actuators:
All actuators are to be of the form with time constant 0.1 seconds,
unless the actual value is available. Appropriate rate and position limits must be specified, and
adhered to.

4. Performance Specifications:
These will be selected by you, unless indicated otherwise.
4.1 Time domain specifications
4.1.1 Step Input for controllers c, d, e

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

%% Part 3: Actuators
% Given initail continuous system (A,B)
A = [ -0.132    0.324    -0.94      0.149  0;...
      -10.614  -1.179     1.0023    0      0; ...
       0.997   -0.00182  -0.259     0      0;...
       0         1        0.34      0      0;...
       0         0        1.0561    0      0]; 
% A B for actuators：make somthing move or operate
B = [0.0069 0.0189; -5.935  1.203; -0.122   -0.614; 0    0; 0    0];
C = diag([1,1,1,1,1]);
D = zeros(5,2); % 5 rwo, 2 col all zero

% Creat new matrix (Why?): Add actuator dynamics
% A & B Matrix with actuator dynamics

% % Remove ψ (last row)
% A(:,5) = [];
% A(5,:) = [];
% B(5,:) = [];

% How to get C&D?
%  s/(s+ 1/tau), where tau = 0.1s
A = [A B; zeros(2,5) [-10 0;0 -10]]
% A_new(:,5) = []
B = [zeros(5,2); [10 0; 0 10]]
C = [diag([1,1,1,1,1,1,1]); zeros(2,7); 0 0 0 0 0 -10 0; 0 0 0 0 0 0 -10]
D = [zeros(7,2); 1 0; 0 1; 10 0; 0 10]
% 7b: Controller gains and the weighting matrices used to calculate them.

% How to get Q & R? Define by us,tsting start from diag([1,1,1,1,1,1])
% Designer selects continuous Q and R, algorithm calculates dQ，dR, M.

Q = diag([25,0,1,50,1,1,1]);
Q = [1e-100        0     0    0       0   0   0; 
         0       1e-100  0    0       0   0   0;
         0         0     1    0       0   0   0;
         0         0     0  1e-100    0   0   0;
         0         0     0    0      30   0   0;
         0         0     0    0       0   3   0;
         0         0     0    0       0   0   5];

% Q = [25 0  0  0  0  0;
%      0  0  0  0  0  0;
%      0  0  1  0  0  0;
%      0  0  0  50 0  0;
%      0  0  0  0  1  0;
%      0  0  0  0  0  1];



R = [1 0; 0 1];

[K, Qd, Rd, M, P, E] = lqrdjv(A, B, Q, R, T)

disp('(7b) Controller Gains K: ')
disp(K)

disp('(7b) Weighting matrix Q: ')
disp(Q)
disp('(7b) Weighting matrix R: ')
disp(R)


%% 7a: Calculate Closed-loop characteristics: eigenvalues, eigenvectors, 
% natural frequencies, damping ratios and time constants as applicable 
% characteristics 
% S-plane eigenvalues
% S-plane is a complex plane with an imaginary and real axis referring to the complex-valued variable  z .
sp_eig = eig(A-B*K);
 [eigenvectors_V, eigenvalues_D] = eig(A-B*K);

disp('Closed-loop characteristics:')
disp('(7a) Eigenvalues (lambda_D) is: ')
disp(eigenvalues_D)
disp('(7a) Eigenvectors (V) is: ')
disp(eigenvectors_V) 

damp(A-B*K)
% disp("S Plane Eigenvalues:");
% disp(sp_eig);

%% Imaginary part of complex number: imag
if imag(sp_eig(1)) == 0
    tao1 = 1/abs(sp_eig(1)); % formula 是什么公式
    ts2_1 = 4/(1/tao1);
    tr_1 = (1 + 1.1*1 + 1.4*1^2) / (1/tao1);
    % hw1 solution
else
    wn1 = abs(sp_eig(1));      % natural frequencies
    wd1 = imag(sp_eig(1));     % damping frequencies
    dr1 = -real(sp_eig(1))/wn1;% damping ratio
    ts2_1 = 4/(dr1*wn1);       % time const
    tr_1 = (1+1.1*dr1+1.4*dr1^2)/(wn1); % ?
end

%%
if imag(sp_eig(2)) == 0
    tao2 = 1/abs(sp_eig(2));
    ts2_2 = 4/(1/tao2); % Settling time 2% Overshoot T_s_2%
    tr_2 = (1+1.1*1+1.4*1^2)/(1/tao2); % Rise Time T_r_90%
else
    wn2 = abs(sp_eig(2));
    wd2 = imag(sp_eig(2));
    dr2 = -real(sp_eig(2))/wn2;
    ts2_2 = 4/(dr2*wn2);
    tr_2 = (1+1.1*dr2+1.4*dr2^2)/(wn2);
end

%%
if imag(sp_eig(3)) == 0
    tao3 = 1/abs(sp_eig(3));
    ts2_3 = 4/(1/tao3);
    tr_3 = (1+1.1*1+1.4*1^2)/(1/tao3);
else
    wn3 = abs(sp_eig(3));
    wd3 = imag(sp_eig(3));
    dr3 = -real(sp_eig(3))/wn3;
    ts2_3 = 4/(dr1*wn3);
    tr_3 = (1+1.1*dr3+1.4*dr3^2)/(wn3);
end

%%
if imag(sp_eig(4)) == 0
    tao4 = 1/abs(sp_eig(4));
    ts2_4 = 4/(1/tao4);
    tr_4 = (1+1.1*1+1.4*1^2)/(1/tao4);
else
    wn4 = abs(sp_eig(4));
    wd4 = imag(sp_eig(4));
    dr4 = -real(sp_eig(4))/wn4;
    ts2_4 = 4/(dr4*wn4);
    tr_4 = (1+1.1*dr4+1.4*dr4^2)/(wn4);
end

%%
if imag(sp_eig(5)) == 0
    tao5 = 1/abs(sp_eig(5));
    ts2_5 = 4/(1/tao5);
    tr_5 = (1+1.1*1+1.4*1^2)/(1/tao5);
else
    wn5 = abs(sp_eig(5));
    wd5 = imag(sp_eig(5));
    dr5 = -real(sp_eig(5))/wn5;
    ts2_5 = 4/(dr5*wn5);
    tr_5 = (1+1.1*dr5+1.4*dr5^2)/(wn5);
end

%%
if imag(sp_eig(6)) == 0
    tao6 = 1/abs(sp_eig(6));
    ts2_6 = 4/(1/tao6);
    tr_6 = (1+1.1*1+1.4*1^2)/(1/tao6);
else
    wn6 = abs(sp_eig(6));
    wd6 = imag(sp_eig(6));
    dr6 = -real(sp_eig(6))/wn6;
    ts2_6 = 4/(dr6*wn6);
    tr_6 = (1+1.1*dr6+1.4*dr6^2)/(wn6);
end


% 4.1 Time domain specifications
x0 = [0; 0; 0; pi/18; 0; 0];
ym = [0; 0; 0; 0; 0; 0];
t_final = 10;
h = 0.01;


%%
% 7(b) Controller gains and the weighting matrices used to calculate them
% Get discrete phi and gamma： 

sysc = ss(A, B, C, D);
sysd = c2d(sysc, h)
phi_cap = sysd.A;
gamma_cap = sysd.B;

%% 7c
% calculate number of frames
frame_number = t_final/h;
time = zeros(frame_number, 1);
data_set1 = zeros(frame_number, 1);
data_set2 = zeros(frame_number, 1);
data_set3 = zeros(frame_number, 1);
data_set4 = zeros(frame_number, 1);
data_set5 = zeros(frame_number, 1);
data_set6 = zeros(frame_number, 1);
data_set7 = zeros(frame_number, 1);
data_set8 = zeros(frame_number, 1);
data_set9 = zeros(frame_number, 1);
data_set10 = zeros(frame_number, 1);
command_set1 = zeros(frame_number, 1);
command_set2 = zeros(frame_number, 1);


% set initial value of x
xk = x0;

% set counter 
counter = 0;

% find u0
uk = K*(ym-x0);

for i=1:frame_number
    % zero order hold
    if(counter == T/h)
        % calculate uk
        uk = [K*(ym - xk)];
        counter = 0;
    end
    % update counter
    counter = counter + 1;

    % perform major computation
    xk1 = phi_cap*xk + gamma_cap*uk;
    yk = C*xk + D*uk;
    
    % store time, command and response (convert to degrees); conversion 
    time(i) = (i-1) * h;
   data1(i) = 180*yk(1, 1)/pi;
   data2(i) = 180*yk(2, 1)/pi;
   data3(i) = 180*yk(3, 1)/pi;
   data4(i) = 180*yk(4, 1)/pi;
   data5(i) = 180*yk(5, 1)/pi;
   data6(i) = 180*yk(6, 1)/pi;
   data7(i) = 180*yk(7, 1)/pi;
   data8(i) = 180*yk(8, 1)/pi;
   data9(i) = 180*yk(9, 1)/pi;
   data10(i) = 180*yk(10, 1)/pi;
    command_set1(i) = 180*uk(1, 1)/pi;
    command_set2(i) = 180*uk(2, 1)/pi;

    % update xk
    xk = xk1;

end

%% plot
figure;
plot(time,data1);
hold on;
plot(time,data2);
hold on;
plot(time,data3);
hold on;
plot(time,data4);
legend("beta", "p", "r", "phi");
title("Closed loop response");
xlabel("time (sec)");
ylabel("states");

figure;
plot(time,data5);
hold on;
plot(time,data6);
hold on;
plot(time,data7);
hold on;
plot(time,data8);
hold on;
plot(time,data9);
hold on;
plot(time,data10);
legend("dac", "drc", "da", "dr", "dadot", "drdot");
title("Closed loop response");
xlabel("time (sec)");
ylabel("actuator states");

figure;
plot(time, command_set1);
hold on;
plot(time, command_set2);
legend("dac", "drc");
title("Closed loop control");
xlabel("time (sec)");
ylabel("control (degrees)");

%% Ask:
% What if by using "impulse(sysc)" and "stepinfo" to solve?
% https://www.mathworks.com/help//ident/ref/lti.impulse.html
% https://www.mathworks.com/help/ident/ref/lti.stepinfo.html#d124e196926
% impulse(sysc)


%%
% my_sys = c2d(sysc, h)
% damp(my_sys)
% 
% my_sys_SS = stepinfo(my_sys,0.02);
% my_sys_SS(2,1)
% my_sys_SS(2,2)

TR_RiseTime = [tr_1 tr_2 tr_3 tr_4 tr_5 tr_6];
fprintf('TR_Rise Time is \n %d %d %d %d %d %d seconds.\r\n',TR_RiseTime)

TS_SettlingTime = [ts2_1  ts2_2  ts2_3  ts2_4  ts2_5  ts2_6];
fprintf('TS_Settling Time is \n %d %d %d %d %d %d seconds.',TS_SettlingTime)











