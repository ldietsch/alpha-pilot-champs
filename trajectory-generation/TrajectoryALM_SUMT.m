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
info.nMPC = 40;
info.Nsteps = 10;
info.substeps = 40; %steps to be used with integrator
info.rho = .1; %multiplier to control inputs in cost function
nMPC = info.nMPC; %N steps ahead in MPC
Ts = 0.01;
info.Ts = Ts; %time step to be used in integrator
info.post_processing = false;
% Q0 = [diag([49 49 0]) zeros(3) zeros(3) zeros(3); zeros(9,3) zeros(9)];
% Q = getQbar(Q0,n,nMPC);
info.R = eye(info.dimM*info.Nsteps);
%% SUMT Initialization

Gates = readtable('Gate Locations.xlsx');
gate.N = Gates{:,1};
gate.x = Gates{:,2};
gate.y = Gates{:,3};
gate.z = Gates{:,4};
gate.current = 0;
gate.last = gate.N(end);

for k = 1:1:gate.last
    if k == 1
        info.x0 = zeros(12,1);
    else
        info.x0 = xref(end-11:end);%previous states
    end
    
    info.xf = [gate.x(k+1),gate.y(k+1),gate.z(k+1)];
    info.motorlimit = 3*2.6978; %approx 3:1 power to weight
    
    
    U0 = ones(m*info.Nsteps,1) * 2.6978;     %create initial input for minimization
    U0 = [U0; Ts];            %add Ts to input vector
    Nsteps = info.Nsteps;
    info.post_processing = true;
    info.two_sample_times = false;
    
    iter = 1;
    info.rp = 1;
    
    
    
    [h, xref] = equalCon(U0, info);
    fnew = costFun(U0, xref, info);
    g = inequalCon(U0, info);
    
    fold = 2 * fnew;
    
    x0 = info.x0;
    xref = predictStates(x0,U0,info);
    info.k_eq = 3;
    
    flam = @(U)costFun(U0, xref, info);
    hlam = @(U)equalCon(U0, info);
    info.lamh = getLamZero(info.x0, U0, flam, hlam, info);
    
    info.lamg = zeros(size(g));
    gam = 3;
    
    options = optimoptions(@fminunc,'Algorithm', 'quasi-newton', ...
        'Display', 'iter'); %,'SpecifyObjectiveGradient', false);
    
    maxiter = 100;
    
    tic
    while abs(fnew-fold) >= 1e-6 && iter<maxiter
        lambda_g = info.lamg;
        lambda_h = info.lamh;
        rp = info.rp;
        
        fold = fnew;  % store last objective function value
        % call fminunc - use "ALM" pseudo-objective function, note that r_p and
        % lambda are passed as parameters, no semi-colon to display results
        [Ustar,phistar,exitflag] = fminunc(@(U)psuedoFun(U, info),U0,options);
        % compute objective and constraints at current xstar
        
        
        [h, xref] = equalCon(Ustar, info);
        f = costFun(Ustar, xref, info);
        g = inequalCon(Ustar, info);
        
        
        fnew = f;
        
        % update lagrange multipliers
        for i = 1:Nsteps+2
            lambda_g(i) = lambda_g(i) + 2 * rp * max( g(i), -lambda_g(i)/(2*rp));
        end
        
        
        for j = 1:3
            lambda_h(j) = lambda_h(j) + 2*rp*h(j);
        end
        iter = iter + 1;     % increment minimization counter
        info.rp = rp * gam; % increase penalty multiplier (gamma = 5, here)
        U0 = Ustar;    % use current xstar as next x0
        info.Ts1 = Ustar(end);
        
    end
    toc/60
    
    x0 = info.x0;
    xref = predictStates(x0,Ustar,info);
    [xi, yi, zi] = extractPos(xref, info.Nsteps+1);
    figure(k)
    plot3(xi,yi,-zi)
end

























