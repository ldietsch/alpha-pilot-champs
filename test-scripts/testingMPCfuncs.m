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
info.dimM = m;
info.Nstates = n;
info.nMPC = 15;
info.Nsteps = 15;
info.substeps = 5; %steps to be used with integrator
info.rho = .1; %multiplier to control inputs in cost function
nMPC = info.nMPC; %N steps ahead in MPC
Ts = 0.01;
info.Ts = Ts; %time step to be used in integrator
info.post_processing = false;
Q0 = [diag([1 1 1]) zeros(3) zeros(3) zeros(3); zeros(9,3) zeros(9)];
Q = getQbar(Q0,n,nMPC);
R = eye(info.dimM*info.nMPC);

x0 = zeros(n,1)*.2;
Uref = [ones(m*(nMPC),1)']';

%% testing main MPC algorithm
info.post_processing = false;
A = []; b = []; Aeq = []; beq = []; 
a_max = 3*g;
lb = zeros(m*nMPC,1); 
ub = (m_tot*a_max/4)*ones(m*nMPC,1);
v_max = 80*1000/60^2; %80 km/h
quad_limits = struct('a_max',a_max,'v_max',v_max);
Umpc = zeros(size(Uref));
options = optimoptions('fmincon','Algorithm','sqp');
j = 1;
k = 1;
l = 1;
x_final_pos = [.5 .5 .5]';
n_total = info.nMPC;
% for i = 1:nMPC
%     info.nMPC = n_total - i + 1;
    nMPC = info.nMPC;
    Q = getQbar(Q0,n,nMPC);
    R = eye(info.dimM*info.nMPC);
    lb = zeros(m*nMPC,1); 
    ub = (m_tot*a_max/4)*ones(m*nMPC,1);
%     if i ~= 1
%         U0 = U(1+m:m*nMPC+m, 1);
%     else
        U0 = [Uref(k:m*nMPC + 1*m - m, 1)];
%     end
    k = k + m;
    j = j + n;
    [Umpc, fval, exitflag, output] = fmincon(@(U)costFunction(U,Q,R,x0,...
    info),U0, A, b, Aeq, beq, lb, ub, ...
    @(U)nlcons(U, x0, x_final_pos, quad_limits, info),options);
%     Umpc(l:l+3) = U(1:m);
    x0 = RK45(x0, Umpc(l:l+3), info);
    l = l + m;
% end


%% Post-processing
info.post_processing = true;
x0 = zeros(n,1)*.2;
results = predictStates(x0, Umpc, info);
[x, y, z] = extractPos(results, info.Nsteps);
plot3(x,y,z,'--')


