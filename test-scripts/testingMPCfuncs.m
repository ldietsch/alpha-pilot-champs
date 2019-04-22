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
info.nMPC = 20;
info.Nsteps = 40;
info.substeps = 5; %steps to be used with integrator
info.rho = .1; %multiplier to control inputs in cost function
nMPC = info.nMPC; %N steps ahead in MPC
Ts = 0.01;
info.Ts = Ts; %time step to be used in integrator
info.post_processing = false;
Q0 = [diag([49 49 0]) zeros(3) zeros(3) zeros(3); zeros(9,3) zeros(9)];
Q = getQbar(Q0,n,nMPC);
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
for i = 1:n_total
    if (nMPC == n_total - i)
       x_final_pos = xf;
    elseif (nMPC > n_total - i) 
           info.nMPC = n_total - i + 1;
           nMPC = info.nMPC;
           Q = getQbar(Q0,n,nMPC);
           R = eye(info.dimM*info.nMPC);
           lb = zeros(m*nMPC,1); 
           ub = (m_tot*a_max/4)*ones(m*nMPC,1);
    end
    %xref will always be the next Nmpc points available
    xr = xref(n + j:n*nMPC + i*n);
    U0 = Uref(k:m*nMPC + i*m - m, 1);
    k = k + m;
    j = j + n;
[U, fval, exitflag, output] = fmincon(@(U)costFunction(U,xr,Q,R,x0,...
    info),U0, A, b, Aeq, beq, lb, ub, ...
    @(U)nlcons(U, x0, x_final_pos, quad_limits, info),options);
    Umpc(l:l+3) = U(1:m);
    x0 = RK4(x0, Umpc(l:l+3),Ts,info);
    l = l + m;
end


%% Post-processing
info.post_processing = true;
x0 = zeros(n,1)*.2;
results = predictStates(x0, Umpc, info);
[x, y, z] = extractPos(results, info.Nsteps);
hold on
plot3(x,y,z,'--')


