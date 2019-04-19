clear
clc
addpath('../discretization');
addpath('../motion-model');
addpath('../post-processing');
addpath('../standard-problem');
addpath('../test-functions');
addpath('../mpc-and-guidance');
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
info.nMPC = 40;
info.Nsteps = 40;
info.substeps = 5; %steps to be used with integrator
info.rho = .1; %multiplier to control inputs in cost function
nMPC = info.nMPC; %N steps ahead in MPC
Ts = 0.01;
info.Ts = Ts; %time step to be used in integrator
info.post_processing = false;
info.two_sample_times = false;
info.Ts1 = Ts;
info.Ts2 = [];
info.getPosition = false;
info.getting_next_state = false;
Q0 = eye(3);
Q = getQbar(Q0,3,nMPC);
R = eye(info.dimM*info.nMPC);
%% testing reference trajectory generation
Nsteps = 40;
x0 = zeros(n,1)*.2;
Uref = [7 0 0 0 zeros(m*(Nsteps-1),1)']';
info.post_processing = true;
xref = predictStates(x0,Uref,info); 
[xi, yi, zi] = extractPos(xref, Nsteps+1);
figure(1)
plot3(xi,yi,zi)


%% testing main MPC algorithm
info.post_processing = false;
A = []; b = []; Aeq = []; beq = []; 
a_max = 3*g;
lb = zeros(m*nMPC,1); 
ub = (m_tot*a_max/4)*ones(m*nMPC,1);
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
for i = 1:n_total-info.nMPC + 1
    if (nMPC == n_total - i + 1)
       x_final_pos = xf;
    elseif (nMPC > n_total - i) 
           info.nMPC = n_total - i + 1;
           nMPC = info.nMPC;
           Q = getQbar(Q0,3,nMPC);
           R = eye(info.dimM*info.nMPC);
           lb = zeros(m*nMPC,1); 
           ub = (m_tot*a_max/4)*ones(m*nMPC,1);
    end
    %xref will always be the next Nmpc points available
    xr = xref(3 + j:3*nMPC + i*3);
    U0 = Uref(k:m*nMPC + i*m - m, 1);
    k = k + m;
    j = j + 3;
[U, fval, exitflag, output] = fmincon(@(U)costFunction(U,U0,xr,Q,R,x0,info),...
    U0, A, b, Aeq, beq, lb, ub, ...
    @(U)nlcons(U, x0, x_final_pos, quad_limits, info),options);
    Umpc(l:l+3) = U(1:m);
    info.getting_next_state = true;
    x0 = predictStates(x0, Umpc(l:l+3),info);
    x0 = x0(13:end);
    info.getting_next_state = false;
    l = l + m;
end


%% Post-processing
output_file = '../input-files/test.xlsx';
info.post_processing = true;
x0 = zeros(n,1)*.2;
results = predictStates(x0, Umpc, info);
[x, y, z] = extractPos(results,info.Nsteps);
pos = extractPosVec(results, info.Nsteps);
results = cell(8,2);
results(:,1) = {'Num_Gates', 'IDs','Nsteps','Substeps','Ts','Xref','Uref','GatePositions'};
results(1,2) = {1};
results(2,2:3) = num2cell(1:2);
results(3,2) = {info.Nsteps};
results(4,2) = {info.substeps};
results(5,2) = {info.Ts};
results(6,2:3*(Nsteps)+1) = num2cell(pos');
results(7,2:m*(Nsteps)+1) = num2cell(Umpc);
results(8,2:4) = num2cell([0.5 0.5 0.5]);
xlswrite(output_file, results);
hold on
plot3(x,y,z,'--')


