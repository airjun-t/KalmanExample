
% https://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling
clear; close all; clc;
M = 5; %cart
m = 1; %pendulum
b = 1;
% I = 0.006;
g = -10;
l = 2;
L = 2;
% p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices

% A = [0      1              0           0;
%      0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
%      0      0              0           1;
%      0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
% B = [     0;
%      (I+m*l^2)/p;
%           0;
%         m*l/p];


A = [0      1               0         0;
     0    -b/M            m*g/M       0;
     0      0               0         1;
     0   -b/(M*L)    (m+M)*g/(M*L)    0];


B = [0;
     1/M;
     0;
     1/(M*L)];

C = [1 0 0 0];

D = [0;
     0];

Q = 0.1 * eye(4); % process (state model) noise covaraince
R = 1; % measurement noise

BF = [B Q 0*B ];

% input to the system is u, disturbance for each state, and measurement
% noise
sysC = ss(A,BF,C,[0 0 0 0 0 R]);



%% Build Kalman Filter
sysFullOutput = ss(A,BF,eye(4),zeros(4,size(BF,2)));

[Kf,P,E] = lqe(A,Q,C,Q,R);

sysKF = ss(A-Kf*C, [B Kf],eye(4), 0*[B Kf]);

%% Setup Sim
dt = 0.01;

T = 50;
t = 0:dt:T;

uD = randn(4,size(t,2));
uN = randn(1,size(t,2));


u = 0*t;
u(100:120) = 100;
u(1500:1520) = -100;

% why is Vd squared
uAUG = [u; Q*uD; R*uN];

figure(1);
[y_sim,t] = lsim(sysC, uAUG,t);
plot(t,y_sim,'b');
hold on;
[x_true,t] = lsim(sysFullOutput, uAUG, t);
[x_est,t] = lsim(sysKF,[u ; y_sim'],t);

plot(t,C*x_est',"k--");
legend('noisy','kalman estimate') 

hold off;
% plot(t,xtrue);



%% Discrete simulation
         
T = 50;
dt = 0.01;
t = 0:dt:T;


F_sys = @(time,state,input) A * state + BF* input;
F_model = @(time,state,input) A * state + B * input;
F_Pdot = @(time,state,input) A * P * A' + Q;


x = zeros(4,length(t));
x_hat = zeros(4,length(t));
y = zeros(1,length(t));
u = zeros(1,length(t));

%Intial Conditions
x0 = [0; 0; 0; 0];
x(:,1) = x0;
x_hat(:,1) = x0;
y(1) = C * x0;
P = zeros(4,4);
u(1) = 0;

for k = 1:length(t)-1                                                                                                           

% to demonstrate basic filter, uncomment lines 133-140ish and comment out 
% and comment out the LQR control lines (146-153)


    if k >=100 && k < 120 
        u(k) = 100;
    elseif k >= 1500 && k <= 1520 
        u(k) = -100;
    else
        u(k) = 0;
    end

    wk = Q*randn(4,1);
    vk = randn(1,1);


    input = [u(k); wk; vk];

    y(k) = C*x(:,k)+ R *vk; % get measurement

    
    % get initial, unfiltered estimate
    x_temp = rk4(F_model, t(k), dt, x_hat(:,k), u(:,k));
    P_temp = rk4(F_Pdot, t(k), dt, P, 0);
    
    [x_next, Kf, P_next] = estimate(P_temp, C, R, x_temp, y(k));


    x_hat(:,k+1) = x_next;
    P = P_next;
  
    
    x(:,k+1) = rk4(F_sys,t(k),dt, x(:,k), input); 
    

end

% 
%% Plot results
axis tight;
figure(1);
plot(t,y,'DisplayName','noisy');
hold on;
plot(t,C*x_hat,'k--','DisplayName','estimate');
% plot(t,C*x,'ro','DisplayName','real');
hold off;
legend();

figure(2);


plot(t,x(2,:),'DisplayName','v true','LineWidth',2);
hold on;

plot(t,x_hat(2,:),'k--','DisplayName','v hat','LineWidth',2);

plot(t,x(3,:),'DisplayName','\theta true','LineWidth',2);
plot(t,x_hat(3,:),'k--','DisplayName','\theta hat','LineWidth',2);

plot(t,x(4,:),'DisplayName','\omega true','LineWidth',2);
plot(t,x_hat(4,:),'k--','DisplayName','\omega hat','LineWidth',2);

title("Unmeasured States");
legend();
hold off;
