%% NUMERICAL EXAMPLE ON CONTRACTIVITY
%  system: one-link robot arm
%  REMARKS
%  Controller design based on Lyapunov's indirect method (First-order Taylor's expansion)

%% Initialization
clear,clc
rng(1);

%% System parameters
global K xstar

n = 4; % dimension of state
m = 1; % dimension of input

J1=0.15; J2=0.2;
F1=0.1; F2=0.15;
mas=0.4; g=9.8; d=0.1; Kc=0.4; Nc=2; % system parameters 

% desired equilibrium point
xstar = [-0.360653451651053 1.185376917919984e-10  1.112599294018070  2.060985191391331e-07]';
ustar = -0.183391898835813; 

%% Data acquisition phase via Simulink
T    = 10;  % number of samples
% Ts   = 0.1; % sampling time 
% Tsim = T*Ts; % duration of simulation
% 
% mag = 0.1;
% x0  = (2*mag).*rand(n,1)-mag + xstar; % initialize the state and control input around the desired equilibrium point to generate data
% 
% sim('data_collection_arm_taylor');
% 
% x  = state.signals.values'; 
% xd = state_deriv.signals.values'; 
% u  = input.signals.values'; 
% 
% X0  = x(:,1:T);
% U0  = u(:,1:T);
% X1  = xd(:,1:T);

load ("Ex1_data2.mat","X0","U0","X1"); % data to generate the controller in the paper and you can generate new code using the previous code

% Data for shifted system
X0tilde = [];
for i=1:T
    X0tilde = [X0tilde X0(:,i)-xstar];
end
U0tilde = U0 - ones(1,T)*ustar;

% Solve the upper bound of Taylor's reminder
delta = 4*abs(cos(X0tilde(1,:))-ones(1,T)) + 2*abs(sin(X0tilde(1,:))-X0tilde(1,:));
Delta = sqrt(T)*[0 0 0 0; 0 max(delta) 0 0; 0 0 0 0; 0 0 0 0];

%% Controller design (First-order Taylor's expansion)
cvx_begin sdp
    variable P1(n,n) symmetric
    variable Y1(T,n)
    variable a 
    variable miu
    P1 >= 10^-6*eye(n);
    a >= 10^-6;
    miu >= 10^-6;
    X0tilde*Y1 == P1;
    [X1*Y1+transpose(X1*Y1)+a*eye(n)+miu*Delta^2*eye(n) transpose(Y1);
    Y1 -miu*eye(T)] <= 0;
cvx_end

G = Y1/P1;
K = U0tilde*G; 

%% Estimate of ROA for shifted closed-loop system
range = 6;
v = -1*range : 0.01 : range;  

[x1 x2] = meshgrid(v);
[e,f] = size(x1);

% Projection of the grid onto a 2D space (\tilde x_1, \tilde x_2) with \tilde x_3 and \tilde x_4 fixed as 0
x3 = zeros(e,f);
x4 = zeros(e,f);

P1inv=inv(P1);

% Vector Q(\tilde x)
Qx1  = cos(x1)-1;
Qx2  = sin(x1);

% Vector \dot \tilde x
M = [0 1 0 0; -Kc/J2 -F2/J2 Kc/J2/Nc 0; 0 0 0 1; -Kc/J1/Nc+K(1)/J1 K(2)/J1 Kc/J1/Nc/Nc+K(3)/J1 -F1/J1+K(4)/J1];% linear part of closed-loop system
N = [0 0; -mas*g*d/J2*cos(xstar(1)) mas*g*d/J2*sin(xstar(1)); 0 0; 0 0];% nonlinear part of closed-loop system
xdot1 = M(1,1)*x1+M(1,2)*x2+M(1,3)*x3+M(1,4)*x4+N(1,1)*Qx1+N(1,2)*Qx2;
xdot2 = M(2,1)*x1+M(2,2)*x2+M(2,3)*x3+M(2,4)*x4+N(2,1)*Qx1+N(2,2)*Qx2;
xdot3 = M(3,1)*x1+M(3,2)*x2+M(3,3)*x3+M(3,4)*x4+N(3,1)*Qx1+N(3,2)*Qx2; 
xdot4 = M(4,1)*x1+M(4,2)*x2+M(4,3)*x3+M(4,4)*x4+N(4,1)*Qx1+N(4,2)*Qx2; 

% H1 = \tilde x'*P1^{-1}*\dot \tilde x
H1 =     P1inv(1,1)*x1.*xdot1 + P1inv(1,2)*x1.*xdot2 + P1inv(1,3)*x1.*xdot3 + P1inv(1,4)*x1.*xdot4...
       + P1inv(2,1)*x2.*xdot1 + P1inv(2,2)*x2.*xdot2 + P1inv(2,3)*x2.*xdot3 + P1inv(2,4)*x2.*xdot4...
       + P1inv(3,1)*x3.*xdot1 + P1inv(3,2)*x3.*xdot2 + P1inv(3,3)*x3.*xdot3 + P1inv(3,4)*x3.*xdot4...
       + P1inv(4,1)*x4.*xdot1 + P1inv(4,2)*x4.*xdot2 + P1inv(4,3)*x4.*xdot3 + P1inv(4,4)*x4.*xdot4;


% H2 = \tilde x'*P1^{-1}*\tilde x
H2 =     P1inv(1,1)*x1.*x1 + P1inv(2,1)*x2.*x1 + P1inv(3,1)*x3.*x1 + P1inv(4,1)*x4.*x1+...
       + P1inv(1,2)*x1.*x2 + P1inv(2,2)*x2.*x2 + P1inv(3,2)*x3.*x2 + P1inv(4,2)*x4.*x2+...
       + P1inv(1,3)*x1.*x3 + P1inv(2,3)*x2.*x3 + P1inv(3,3)*x3.*x3 + P1inv(4,3)*x4.*x3+...
       + P1inv(1,4)*x1.*x4 + P1inv(2,4)*x2.*x4 + P1inv(3,4)*x3.*x4 + P1inv(4,4)*x4.*x4;

H3 = zeros(e,f); % set of the points for which the Lyapunov derivative is positive inside the Lyapunov sublevel set H2
index   = []; % record the number of points in H3

gamma = 16.3; % range of Lyapunov sublevel set H2; modify gamma until H3 is empty

for i= 1:e
    for j= 1:f
        if H2(i,j) <= gamma
            if H1(i,j) > 0
                index = [index [i j]'];
            end
        else
            continue;
        end
    end
end

number_index = size(index);
for i=1 : number_index(2)
    H3(index(1,i),index(2,i)) = 1;
end

cond1 = H1 < 0; % Lyapunov derivative is negative
cond1 = double(cond1);  
cond1(cond1 == 0) = NaN;
colormap([0.5 0.5 0.5;0 0 1;1 0 0])% gray, blue and red  
surf(x1,x2,cond1,ones(length(v))); % gray
shading interp

hold on
cond2 = H2 <= gamma; % Lyapunov sublevel set
cond2 = double(cond2);
cond2(cond2 == 0) = NaN;
surf(x1,x2,cond2,ones(length(v))+1); % blue
shading interp

hold on
cond3 = H3 > 0; 
cond3 = double(cond3);
cond3(cond3 == 0) = NaN;
surf(x1,x2,cond3,ones(length(v))+2); % red
shading interp

axis([-range range -range range])
view(0,90) 
set(gca,'xtick',-range:2:range) 
set(gca,'ytick',-range:2:range) 

%% Evaluation of obtained controller via ODE45 function
mx0 = 1;   % magnitude of initial conditions
x0  = (2*mx0).*rand(n,1)-mx0;
tspan = [0,100];  % duration of the simulation

[t,x] = ode45(@arm,tspan,x0);
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

function dxdt = arm(t,x)  % shifted system
    global K xstar 
    J1=0.15; J2=0.2;
    F1=0.1; F2=0.15;
    mas=0.4; g=9.8; d=0.1; Kc=0.4; Nc=2;

    dxdt = zeros(4,1);
    dxdt(1) = x(2);
    dxdt(2) = -Kc/J2*x(1) - F2/J2*x(2) + Kc/J2/Nc*x(3) - mas*g*d/J2*cos(xstar(1))*(cos(x(1))-1) + mas*g*d/J2*sin(xstar(1))*sin(x(1));
    dxdt(3) = x(4);
    dxdt(4) = -Kc/J1/Nc*x(1) + Kc/J1/Nc/Nc*x(3) - F1/J1*x(4) + 1/J1*K*x;
end