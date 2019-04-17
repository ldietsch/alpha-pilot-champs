clear all
close all
clc

%% adding necessary directories
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
info.nMPC = 20;
info.substeps = 5; %steps to be used with integrator
info.rho = .1; %multiplier to control inputs in cost function
nMPC = info.nMPC; %N steps ahead in MPC
info.post_processing = false;

%% testing single-shooting nmpc function
% guidance trajectory goes here
[uMPC,~,exitflag,~] = ssnmpc(input_file,info);


%% Post-processing
info.post_processing = true;
x0 = zeros(n,1)*.2;
results = predictStates(x0, uMPC, info);
[x, y, z] = extractPos(results, info.Nsteps);
hold on
plot3(x,y,z,'--')


