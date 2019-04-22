function results = funcTrajectoryALM_SUMT_FixedTime(input_file, runNum)
addpath('../discretization');
addpath('../motion-model');
addpath('../post-processing');
addpath('../standard-problem');
addpath('../test-functions');

%% defining problem parameters
Mbig = 1;
Msmall = 0.1;
m_tot = Mbig + Msmall;
Radius = 0.2;
l = 0.1;
J = getMoments(Mbig, Msmall, Radius, l);
k_tau = 1;
g = 9.81;
tol = 1e-6;
info = struct('m',m_tot,'g',g,'J',J,'l',l,'k_tau',k_tau,'tol',tol);
%% reading in data from input file
[~,~,params] = xlsread(input_file,'params');
info.Nsteps = params{runNum+1,2};
info.substeps = params{runNum+1,4};
info.Ts = params{runNum+1,3};
n = 12; %Number of states
m = 4; %Number of inputs
info.dimM = m;
info.Nstates = n;
info.rho = .1; %multiplier to control inputs in cost function
info.post_processing = false;
info.vel_upper = 16.66/sqrt(3);
info.R = eye(info.dimM*info.Nsteps);


%% SUMT Initialization

Gates = readtable('Gate Locations.xlsx');
% Gates = readtable('Gate Locations_debug.xlsx');
gate.N = Gates{:,1};
gate.x = Gates{:,2};
gate.y = Gates{:,3};
gate.z = Gates{:,4};
gate.current = 0;
gate.last = 5; % note it's fixed to 5

Gate_ID = [];
Ts_total = [];
X_ref_total = [];
Ustar_total = [];
Gatepos_total = [];
gate_x = []; gate_y = []; gate_z = [];
k = 1;
while k <= gate.last %this should be read as the for-loop. Had to adjust to
    %make it work with step adjustment
    if k == 1 
        info.x0 = zeros(12,1);     
    else
        info.x0 = xref(end-11:end);%previous states
    end
    
    nextGate = [gate.x(k+1),gate.y(k+1),gate.z(k+1)];
    info.xf = nextGate;
        info.psif = [];
        info.k_eq = 3;

    info.motorlimit = 3*2.6978; %approx 3:1 power to weight
    
    gate_x = [gate_x, gate.x(k+1)]; gate_y = [gate_y, gate.y(k+1) ]; gate_z = [gate_z, gate.z(k+1) ];
    
    
    U0 = ones(m*info.Nsteps,1) * 2.6978;     %create initial input for minimization

    info.post_processing = true;
    info.two_sample_times = false;
    
    iter = 1;
    info.rp = 1;
    
    
    
    [h, xref] = equalConFixedTime(U0, info);
    fnew = costFunFixedTime(U0, xref, info);
    g = inequalConFixedTime(U0, info);
    
    fold = 2 * fnew;
    
    x0 = info.x0;
    xref = predictStates(x0,U0,info);

    
    flam = @(U)costFunFixedTime(U0, xref, info);
    hlam = @(U)equalConFixedTime(U0, info);
    info.lamh = getLamZero(info.x0, U0, flam, hlam, info);
    info.lamg = zeros(size(g));
    lambda_g = info.lamg;
    lambda_h = info.lamh;

    gam = 3;
    
    options = optimoptions(@fminunc,'Algorithm', 'quasi-newton', ...
        'Display', 'iter'); %,'SpecifyObjectiveGradient', false);
    
    maxiter = 100;
    
    tic
    while abs(fnew-fold) >= tol && iter<maxiter
        rp = info.rp;
        
        fold = fnew;  % store last objective function value
        % call fminunc - use "ALM" pseudo-objective function, note that r_p and
        % lambda are passed as parameters, no semi-colon to display results
        [Ustar,phistar,exitflag] = fminunc(@(U)psuedoFunFixedTime(U, info),...
            U0,options);
        % compute objective and constraints at current xstar
        
        
        [h, xref] = equalConFixedTime(Ustar, info);
        f = costFunFixedTime(Ustar, xref, info);
        g = inequalConFixedTime(Ustar, info);
        
        
        fnew = f;
        
  
            lambda_g = lambda_g + 2 * rp * max( g, -lambda_g/(2*rp));

            lambda_h = lambda_h + 2.*rp.*h;

        iter = iter + 1;     % increment minimization counter
        info.rp = rp * gam; % increase penalty multiplier (gamma = 5, here)
        U0 = Ustar;    % use current xstar as next x0

        
    end

        x0 = info.x0;
        xref = predictStates(x0,Ustar,info);
        Gate_ID = [Gate_ID, k];
        X_ref_total = [X_ref_total; xref];
        Ustar_total = [Ustar_total; Ustar];
        Gatepos_total = [Gatepos_total, info.xf];
        toc/60
        redo_traj = false;
        k = k + 1;
    
    
    
end

results = num2cell([runNum Ustar_total']);

end


























