clear
clc
addpath('../discretization');
addpath('../motion-model');
addpath('../post-processing');
addpath('../standard-problem');
addpath('../test-functions');
%% defining problem parameters
Mbig = 1;
Msmall = 0.1;
m_tot = Mbig + Msmall;
R = 0.2;
l = 0.1;
J = getMoments(Mbig, Msmall, R, l);
k_tau = 1;
g = 9.81;
tol = 1e-4;
info = struct('m',m_tot,'g',g,'J',J,'l',l,'k_tau',k_tau,'tol',tol);
n = 12; %Number of states
m = 4; %Number of inputs
Nmpc = 5; %N steps ahead in MPC
Ts = 0.01;
Q0 = [eye(3) zeros(3) zeros(3) zeros(3); zeros(9,3) zeros(9)];
Q = getQbar(Q0,n,Nmpc);

%% testing reference trajectory generation
Nsteps = 40;
x0 = zeros(n,1)*.2;
Uref = [0 0 0 0 zeros(m*(Nsteps-1),1)']';
xref = predictStates(x0,Uref,Ts,Nsteps,info); 
position = extractPos(xref, Nsteps);
figure(1)
plot3(position(1,:),position(2,:),position(3,:))
%% testing main MPC algorithm
A = []; b = []; Aeq = []; beq = []; 
a_max = 3*g;
lb = zeros(m*Nmpc,1); 
ub = (m_tot*a_max/4)*ones(m*Nmpc,1);
v_max = 80*1000/60^2; %80 km/h
quad_limits = struct('a_max',a_max,'v_max',v_max);
n_total = Nsteps;
Umpc = zeros(size(Uref));
options = optimoptions('fmincon','Algorithm','sqp');
j = 1;
k = 1;
l = 1;
x_final_pos = [];
xf = xref(Nsteps*12 + 1: Nsteps*12 + 3);
for i = 1:n_total
    if (Nmpc == n_total - i)
       x_final_pos = xf;
    elseif (Nmpc > n_total - i) 
           Nmpc = n_total - i + 1;
           Q = getQbar(Q0,n,Nmpc);
           lb = zeros(m*Nmpc,1); 
           ub = (m_tot*a_max/4)*ones(m*Nmpc,1);
    end
    %xref will always be the next Nmpc points available
    xr = xref(n + j:n*Nmpc + i*n);
    U0 = Uref(k:m*Nmpc + i*m - m, 1);
    k = k + m;
    j = j + n;
[U, fval, exitflag, output] = fmincon(@(U)costFunction(U,xr,Q,x0,...
    Ts,Nmpc,info),U0, A, b, Aeq, beq, lb, ub, ...
    @(U)nlcons(U, x0, x_final_pos, Ts, Nmpc, quad_limits, info),options);
    Umpc(l:l+3) = U(1:m);
    x0 = RK4(x0, Umpc(l:l+3),Ts,info);
    l = l + m;
end


%% Post-processing
steps = length(Umpc)/m;
x0 = zeros(n,1)*.2;
results = predictStates(x0,Umpc, Ts, steps, info);
[x, y, z] = extractPos(results, steps);
figure(2)
plot3(x,y,z)


