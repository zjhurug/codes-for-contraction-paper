%% NUMERICAL EXAMPLE ON CONTRACTIVITY
%  system: one-link robot arm
%  REMARKS
%  In this example we pick Q(x)=cox(x1) or Q(x)=[cox(x1) x1^2 sin(x2)]'

%% Initialization
clear,clc
rng(1);

%% System parameters
global K

n = 4; % dimension of state
m = 1; % dimension of input

J1=0.15; J2=0.2;
F1=0.1; F2=0.15;
mas=0.4; g=9.8; d=0.1; Kc=0.4; Nc=2; % system parameters 

%% Controllability test
Abar = [0 1 0 0; -Kc/J2 -F2/J2 Kc/J2/Nc 0; 0 0 0 1; -Kc/J1/Nc 0 Kc/J1/Nc/Nc -F1/J1];
B = [0 0 0 1/J1]';
rank(ctrb(Abar, B))

%% Data acquisition phase via Simulink
T    = 10;  % number of samples
% Ts   = 0.1; % sampling interval
% Tsim = T*Ts; % duration of simulation
% 
% mag = 0.1; % magnitude of initial conditions
% x0  = (2*mag).*rand(n,1)-mag; % initial state
% 
% sim('data_collection_arm');
% 
% x  = state.signals.values'; 
% xd = state_deriv.signals.values'; 
% u  = input.signals.values'; 
% 
% X0  = x(:,1:T);
% U0  = u(:,1:T);
% X1  = xd(:,1:T);

load("Ex1_data1.mat","X0","U0","X1"); % data to generate the controller in the paper and you can generate new code using the previous code

s = 5; % dimension of Z(x)
Z0  = [X0;cos(X0(1,:))];
RQ = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % Bound for Jacobian of Q(x)

% s = 7; % dimension of Z(x)
% Z0  = [X0;cos(X0(1,:));X0(1,:).^2;sin(X0(2,:))]; % more nonlinearies in Q(x)
% w = 1;
% RQ = [sqrt(4*w^2+1) 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0]; % Bound for Jacobian of Q(x)

rank([U0; Z0]) % full row rank check

%% Controller design (CONTRACTIVITY)
r = n; % column number of RQ

cvx_begin sdp
    variable P1(n,n) symmetric
    variable Y1(T,n)
    variable G2(T,s-n)
    variable a 
    P1 >= 0*eye(n);
    a >= 0;
    Z0*Y1 == [P1;zeros(s-n,n)];
    [X1*Y1+transpose(X1*Y1)+a*eye(n) X1*G2 P1*RQ;
        transpose(X1*G2) -eye(s-n) zeros(s-n,r);
        transpose(P1*RQ) zeros(r,s-n) -eye(r)] <= 0;
    Z0*G2 == [zeros(n,s-n);eye(s-n)];
cvx_end

G1 = Y1/P1;
G  = [G1 G2];
K  = U0*G; 
closed = X1*G;

%% Evaluation of obtained controller via ODE45 function
mx0 = 1;   % magnitude of initial conditions
x0  = (2*mx0).*rand(n,1)-mx0;
tspan = [0,100];  % duration of simulation

[t,x] = ode45(@arm,tspan,x0);

x_star = x(end,:)';
u_star = K*[x_star;cos(x_star(1))]; % equlibrium point of closed-loop system
% u_star = K*[x_star;cos(x_star(1));x_star(1)^2;sin(x_star(2))]; % equlibrium point of closed-loop system

figure
plot(t,x(:,1),'r','LineWidth',1);
hold on;
plot(t,x(:,2),'b','LineWidth',1);
hold on;
plot(t,x(:,3),'g','LineWidth',1);
hold on;
plot(t,x(:,4),'k','LineWidth',1);
xlabel('t');
legend('x(1) ','x(2)','x(3) ','x(4)');

function dxdt = arm(t,x)  
    global K
    J1=0.15; J2=0.2;
    F1=0.1; F2=0.15;
    mas=0.4; g=9.8; d=0.1; Kc=0.4; Nc=2;

    dxdt = zeros(4,1);
    dxdt(1) = x(2);
    dxdt(2) = -Kc/J2*x(1) - F2/J2*x(2) + Kc/J2/Nc*x(3) - mas*g*d/J2*cos(x(1));
    dxdt(3) = x(4);
    dxdt(4) = -Kc/J1/Nc*x(1) + Kc/J1/Nc/Nc*x(3) - F1/J1*x(4) + 1/J1*K*[x;cos(x(1))];
    % dxdt(4) = -Kc/J1/Nc*x(1) + Kc/J1/Nc/Nc*x(3) - F1/J1*x(4) + 1/J1*K*[x;cos(x(1));x(1)^2;sin(x(2))]; % more nonlinearies in Q(x)
end