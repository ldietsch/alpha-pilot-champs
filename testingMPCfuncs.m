clear
clc
Mbig = 1;
Msmall = 0.1;
m_tot = Mbig + Msmall;
R = 0.2;
l = 0.1;
J = getMoments(Mbig, Msmall, R, l);
k_tau = 1;
g = 9.81;
info = struct('m',m_tot,'g',g,'J',J,'l',l,'k_tau',k_tau);
n = 12; %Number of states
N = 20; %N steps ahead in MPC
Ts = 0.05;
Q0 = [eye(3) zeros(3) zeros(3) zeros(3); zeros(9,3) zeros(9)];
Q = getQbar(Q0,n,N);

%% testing the propagation of states
m = 4;
U = ones(m*N,1);
x0 = ones(n,1)*.2;
X = predictStates(x0, U, Ts, N, info);

%% testing reference trajectory generation
xf = [10;10;10];
int = 4;
xref = genRefTraj(x0(1:3),xf,int,N);
plot3(xref(:,1),xref(:,2),xref(:,3))

%% testing main MPC algorithm

A = []; b = []; Aeq = []; beq = [];
options = optimoptions('fmincon','Algorithm','sqp');
[U, fval, exitflag, output] = fmincon(@(U)costFunction(U,Xref,Q,X0,Ts,N),...
    U0,A,b,Aeq,beq,lb,ub,@(U)nlcons(U, X0, x_final, N),options);


