%% NUNERICAL EXAMPLE ON CONTRACTIVITY
% system: bilinear chemical reactor
%  REMARKS
%  In this example we pick Q(\xi)=x1*u=\xi(1)*\xi(3).

%% Initialization
clear,clc
rng(1);

%% System parameters
global K

n = 2; % dimension of state
m = 1; % dimension of input
R = 4; % dimension of Z(\xi)

%% Data acquaition phase
T    = 10;  % number of samples
Ts   = 0.1; % sampling interval for data acquition 
Tsim = T*Ts; % simulation time

mag = 0.1; % magnitude input and initial conditions
x0  = (2*mag).*rand(n+m,1)-mag; % initial state

sim('data_collection_chemical_reactor');

x  = state.signals.values'; 
xd = state_deriv.signals.values'; 
v  = input.signals.values'; 

X0  = x(:,1:T);
V0  = v(:,1:T);
Z0  = [X0; X0(1,:).*X0(3,:)];
X1  = [xd(:,1:T)];

rank([V0; Z0]) % full row rank check

%% Controller design (CONTRACTIVITY)
RQ = [0.1 0 0; 0 0 0; 0 0 0.1]; % Bound for R_Q
r = n+m; % column number of RQ, i.e., dimension of extended state \xi

cvx_begin sdp
    variable P1(n+m,n+m) symmetric
    variable Y1(T,n+m)
    variable G2(T,R-n-m)
    variable a 
    P1 >= 0*eye(n+m);
    a >= 0;
    Z0*Y1 == [P1;zeros(R-n-m,n+m)];
    [X1*Y1+transpose(X1*Y1)+a*eye(n+m) X1*G2 P1*RQ;
        transpose(X1*G2) -eye(R-n-m) zeros(R-n-m,r);
        transpose(P1*RQ) zeros(r,R-n-m) -eye(r)] <= 0;
    Z0*G2 == [zeros(n+m,R-n-m);eye(R-n-m)];
cvx_end

G1 = Y1/P1;
G  = [G1 G2];
K  = V0*G; 
closed = X1*G;

%% Estimate of ROA for closed-loop system
range = 4;
v = -range:0.001:range;  
[y1 y2] = meshgrid(v);
[e,f] = size(y1);

% Projection of the grid onto a 2D space (\xi_1,\xi_3) with \xi_2 fixed as 0
x1=y1;
x3=y2;
x2=zeros(e,f);

P1inv=inv(P1);

% Vector Q(\xi)
Qx = x1.*x3;

% Vector \dot \xi
M = X1*G1;
N = X1*G2;
xidot1 = M(1,1)*x1+M(1,2)*x2+M(1,3)*x3+N(1,1)*Qx; 
xidot2 = M(2,1)*x1+M(2,2)*x2+M(2,3)*x3+N(2,1)*Qx; 
xidot3 = M(3,1)*x1+M(3,2)*x2+M(3,3)*x3+N(3,1)*Qx; 

% H1 = \dot \xi'*P1^{-1}*\xi 
H1 =     P1inv(1,1)*x1.*xidot1 + P1inv(1,2)*x1.*xidot2 + P1inv(1,3)*x1.*xidot3...
       + P1inv(2,1)*x2.*xidot1 + P1inv(2,2)*x2.*xidot2 + P1inv(2,3)*x2.*xidot3...
       + P1inv(3,1)*x3.*xidot1 + P1inv(3,2)*x3.*xidot2 + P1inv(3,3)*x3.*xidot3;

% H2 = \xi'*P1^{-1}*\xi
H2 =     P1inv(1,1)*x1.*x1 + P1inv(2,1)*x2.*x1 + P1inv(3,1)*x3.*x1...
       + P1inv(1,2)*x1.*x2 + P1inv(2,2)*x2.*x2 + P1inv(3,2)*x3.*x2...
       + P1inv(1,3)*x1.*x3 + P1inv(2,3)*x2.*x3 + P1inv(3,3)*x3.*x3;

H3 = zeros(e,f); % set of the points for which the Lyapunov derivative is positive inside the Lyapunov sublevel set H2
index   = []; % record the number of points in H3

gamma = 0.06; % range of Lyapunov sublevel set H2; modify gamma until H3 is empty

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
surf(x1,x3,cond1,ones(length(v))); % gray
shading interp

hold on
cond2 = H2 <= gamma; % Lyapunov sublevel set
cond2 = double(cond2);
cond2(cond2 == 0) = NaN;
surf(x1,x3,cond2,ones(length(v))+1); % blue
shading interp

hold on
cond3 = H3 > 0; 
cond3 = double(cond3);
cond3(cond3 == 0) = NaN;
surf(x1,x3,cond3,ones(length(v))+2); % red
shading interp

axis([-0.4 0.4 -range range])
view(0,90) 
xlabel('\xi_1') 
ylabel('\xi_3') 
set(gca,'xtick',-range:0.2:range) 
set(gca,'ytick',-range:2:range) 

%% Evaluation of obtained controller via ODE45 function
mx0 = 0.1;   % magnitude of initial conditions
x0  = (2*mx0).*rand(n+m,1)-mx0;
tspan = [0,10];  % duration of the simulation

[t,x] = ode45(@reactor,tspan,x0);
figure
plot(t,x(:,1),'r');
hold on;
plot(t,x(:,2),'b');
hold on;
plot(t,x(:,3),'k');
legend('x(1) ','x(2)','u');

function dxdt = reactor(t,x)
    global K
    dxdt = zeros(3,1);
    dxdt(1) = 4.25*x(1) + x(2) - 0.25*x(3) - x(1)*x(3);
    dxdt(2) = -6.25*x(1) - 2*x(2);
    dxdt(3) = K*[x; x(1)*x(3)];
end