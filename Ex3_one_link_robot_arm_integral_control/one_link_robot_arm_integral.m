%% NUMERICAL EXAMPLE ON CONTRACTIVITY
% system: one-link robot arm
% control target: constant output tracking via integral control
%                 constant disturbance: [0.1 0.2 0.3 0.4]
%                 constant reference: r=pi/3  
%  REMARKS
%  In this example we pick Q(x)=cox(x1) or Q(x)=[cox(x1) sin(x2)]'   

%% Initialization
clear,clc
rng(1);

%% System parameters
global K

n = 4; % dimension of state
m = 1; % dimension of input
p = 1; % dimension of output

J1=0.15; J2=0.2;
F1=0.1; F2=0.15;
mas=0.4; g=9.8; d=0.1; Kc=0.4; Nc=2; % system parameters

track_r = pi/3; % track reference

noise = [0.1 0.2 0.3 0.4]; % constant disturbance

%% Data acquisition phase via Simulink
T    = 10;  % number of samples
Ts   = 0.1; % sampling interval
Tsim = T*Ts; % duration of simulation

mag = 0.1; % magnitude of initial conditions
x0  = (2*mag).*rand(n+p,1)-mag; % initial state

sim('data_collection_arm_integral');

x  = state.signals.values'; 
xd = state_deriv.signals.values'; 
u  = input.signals.values'; 

X0  = x(:,1:T);
Z1  = xd(:,1:T);
U0  = u(:,1:T);
M = ones(1,T);

s = 5; % dimension of Z(x)
Z0  = [X0;cos(X0(1,:))];
RQ = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % Bound for Jacobian of Q(x)

% s = 6; % dimension of Z(x)
% Z0  = [X0;cos(X0(1,:));sin(X0(2,:))]; % more nonlinearies in Q(x)
% RQ = [1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % Bound for Jacobian of Q(x)

rank([U0; Z0; M]) % full row rank check

%% Controller design (CONTRACTIVITY)
r = n; % column number of RQ

cvx_begin sdp
    variable P1(n+p,n+p) symmetric
    variable Y1(T,n+p)
    variable G2(T,s-n)
    variable a 
    P1 >= 0*eye(n+p);
    a >= 0;
    Z0*Y1 == [P1;zeros(s-n,n+p)];
    [Z1*Y1+transpose(Z1*Y1)+a*eye(n+p) Z1*G2 P1*RQ;
        transpose(Z1*G2) -eye(s-n) zeros(s-n,r);
        transpose(P1*RQ) zeros(r,s-n) -eye(r)] <= 0;
    Z0*G2 == [zeros(n+p,s-n);eye(s-n)];
    M*[Y1 G2] == zeros(1,s+p);
cvx_end

G1 = Y1/P1;
G  = [G1 G2];
K  = U0*G; 
closed = Z1*G;

%% Evaluation of obtained controller via ODE45 function
mx0 = 1;   % magnitude of initial conditions
x0  = (2*mx0).*rand(5,1)-mx0;
tspan = [0,60];  % duration of the simulation
r = pi/3;

[t,x] = ode45(@arm_integral,tspan,[x0; r]);
figure
plot(t,x(:,1),'r','LineWidth',1);
hold on;
plot(t,x(:,2),'b','LineWidth',1);
hold on;
plot(t,x(:,3),'g','LineWidth',1);
hold on;
plot(t,x(:,4),'k','LineWidth',1);
hold on;
plot(t,x(:,6),'--r','LineWidth',1); % plot track reference r
xlabel('t');
ylabel('x');
legend('x_1 ','x_2','x_3','x_4','r');

function dxdt = arm_integral(t,x)  
    global K
    J1=0.15; J2=0.2;
    F1=0.1; F2=0.15;
    mas=0.4; g=9.8; d=0.1; Kc=0.4; Nc=2;
    noise=[0.1 0.2 0.3 0.4];
    track_r = pi/3;

    dxdt = zeros(6,1);
    dxdt(1) = x(2) + noise(1);
    dxdt(2) = -Kc/J2*x(1) - F2/J2*x(2) + Kc/J2/Nc*x(3) - mas*g*d/J2*cos(x(1)) + noise(2);
    dxdt(3) = x(4) + noise(3);
    dxdt(4) = -Kc/J1/Nc*x(1) + Kc/J1/Nc/Nc*x(3) - F1/J1*x(4) + 1/J1*K*[x(1:5);cos(x(1))]  + noise(4);
    % dxdt(4) = -Kc/J1/Nc*x(1) + Kc/J1/Nc/Nc*x(3) - F1/J1*x(4) + 1/J1*K*[x(1:5);cos(x(1));sin(x(2))]  + noise(4);% more nonlinearies in Q(x)
    dxdt(5) = x(1) - track_r;
    dxdt(6) = 0;
end