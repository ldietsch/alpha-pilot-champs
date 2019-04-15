clear, clc
% testing getting lam0 for ALM in SUMT
%% adding all folders necessary for compilation
addpath('../discretization');
addpath('../motion-model');
addpath('../post-processing');
addpath('../standard-problem');
addpath('../test-functions');
addpath('../SUMT-example');

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
m = 4; %Number of decision variables
info.dimM = m;
info.Nstates = n;
info.nMPC = 15;
info.Nsteps = 20; %in reality, will be a lot bigger
info.substeps = 5;
info.rho = .1; %multiplier to control inputs in cost function
nMPC = info.nMPC; %N steps ahead in MPC
info.optim_sol = true;
Ts = 0.01;
info.Ts = Ts; %time step to be used in integrator
info.post_processing = false;
info.k_eq = 3; % the number of equality constraints
Q0 = [diag([1 1 1]) zeros(3) zeros(3) zeros(3); zeros(9,3) zeros(9)];
Q = getQbar(Q0,n,nMPC);
R = eye(info.dimM*info.Nsteps);
x_final = [10; 10; 10];
x0 = zeros(n,1);
U0 = [1; zeros(info.dimM*info.Nsteps-1,1); info.Ts];

f = @(U)costFunction(U, Q, R, x0, info);
h = @(U)equalityCons(U, x0, x_final, info);

lam0 = getLamZero(x0,U0,f,h,info)