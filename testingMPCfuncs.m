clear
clc
%% defining problem parameters
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
Nmpc = 10; %N steps ahead in MPC
Ts = 0.1;
Q0 = [eye(3) zeros(3) zeros(3) zeros(3); zeros(9,3) zeros(9)];
Q = getQbar(Q0,n,Nmpc);

%% testing the propagation of states
m = 4;
U = ones(m*Nmpc,1);
x0 = ones(n,1)*.2;
X = predictStates(x0, U, Ts, Nmpc, info);

%% testing reference trajectory generation
xf = [10;10;10];
int = 4;
xref = genRefTraj(x0(1:3),xf,int,Nmpc);
plot3(xref(:,1),xref(:,2),xref(:,3))

%% testing main MPC algorithm
Xref = includeAll(xref);
A = []; b = []; Aeq = []; beq = []; lb = []; ub = [];
a_max = 27;
v_max = 80*1000/60^2; %60 km/h
quad_limits = struct('a_max',a_max,'v_max',v_max);
n_total = Nmpc*int;
Uref = ones(m*n_total,1);
Umpc = zeros(size(Uref));
options = optimoptions('fmincon','Algorithm','sqp');
j = 1;
k = 1;
l = 1;
x_final_pos = [];
x0 = [.2 .2 .2 0 0 0 0 0 0 0 0 0]';
for i = 1:n_total
    if (Nmpc == n_total - i)
       x_final_pos = xf;
    elseif (Nmpc > n_total - i) 
           Nmpc = n_total - i+1;
           Q = getQbar(Q0,n,Nmpc);
    end
    %xref will always be the next Nmpc points available
    xref = Xref(j:n*Nmpc+i*n);
    %exclude the first point x_k
    xref = xref(13:end);
    U0 = Uref(k:m*Nmpc+i*m-m, 1);
    k = k + m;
    j = j + n;
[U, fval, exitflag, output] = fmincon(@(U)costFunction(U,xref,Q,x0,Ts,Nmpc,info),...
    U0, A, b, Aeq,beq,lb,ub, ...
    @(U)nlcons(U, x0, x_final_pos, Ts, Nmpc, quad_limits, info),options);
    Umpc(l:l+3) = U(1:m);
    x0 = forwardEuler(x0, Umpc(l:l+3),Ts,info);
    l = l + m;
end

