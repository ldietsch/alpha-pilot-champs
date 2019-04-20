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
tol = 1e-6;
info = struct('m',m_tot,'g',g,'J',J,'l',l,'k_tau',k_tau,'tol',tol);
n = 12; %Number of states
m = 4; %Number of inputs
info.dimM = m;
info.Nstates = n;
info.nMPC = 5;
info.Nsteps = 5;
info.substeps = 40; %steps to be used with integrator
info.rho = 1e-3; %multiplier to control inputs in cost function
nMPC = info.nMPC; %N steps ahead in MPC
Ts = 0.01;
info.Ts = Ts; %time step to be used in integrator
info.post_processing = false;
info.vel_upper = 16.66/sqrt(3);
% Q0 = [diag([49 49 0]) zeros(3) zeros(3) zeros(3); zeros(9,3) zeros(9)];
% Q = getQbar(Q0,n,nMPC);
info.R = eye(info.dimM*info.Nsteps);
%% SUMT Initialization

% Gates = readtable('Gate Locations.xlsx');
Gates = readtable('Gate Locations_debug.xlsx');
gate.N = Gates{:,1};
gate.x = Gates{:,2};
gate.y = Gates{:,3};
gate.z = Gates{:,4};
gate.current = 0;
gate.last = gate.N(end);

Gate_ID = [];
Ts_total = [];
X_ref_total = [];
Ustar_total = [];
Gatepos_total = [];
gate_x = []; gate_y = []; gate_z = [];
redo_traj = false; %initialize boolean
k = 1;
while k <= gate.last %this should be read as the for-loop. Had to adjust to
    %make it work with step adjustment
    if k == 1 && ~redo_traj
%         info.x0 = zeros(12,1);
          info.x0 = [gate.x(k);gate.y(k);gate.z(k); 10; 0; 1;
          0;0;0;0;0;0]; %for debugging only
    elseif redo_traj
        xref = predictStates(x0,Ustar,info);
        info.R = eye(info.dimM*info.Nsteps);        
    else
        info.x0 = xref(end-11:end);%previous states
    end
    
    nextGate = [gate.x(k+1),gate.y(k+1),gate.z(k+1)];
    info.xf = nextGate;
%     if k + 2 < gate.last
%         afterNextGate = [gate.x(k+2),gate.y(k+2),gate.z(k+2)];
%         info.psif = findFinalHeading(nextGate, afterNextGate); 
%         info.k_eq = 4;
%     else 
        info.psif = [];
        info.k_eq = 3;
%     end
    info.motorlimit = 3*2.6978; %approx 3:1 power to weight
    
    gate_x = [gate_x, gate.x(k+1)]; gate_y = [gate_y, gate.y(k+1) ]; gate_z = [gate_z, gate.z(k+1) ];
    
    
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

    
    flam = @(U)costFun(U0, xref, info);
    hlam = @(U)equalCon(U0, info);
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
        [Ustar,phistar,exitflag] = fminunc(@(U)psuedoFun(U, info),U0,options);
        % compute objective and constraints at current xstar
        
        
        [h, xref] = equalCon(Ustar, info);
        f = costFun(Ustar, xref, info);
        g = inequalCon(Ustar, info);
        
        
        fnew = f;
        
        % update lagrange multipliers
%         for i = 1:8*Nsteps+2 % don't actually need recursion
            lambda_g = lambda_g + 2 * rp * max( g, -lambda_g/(2*rp));
%         end
        
        
%         for j = 1:3
            lambda_h = lambda_h + 2.*rp.*h;
%         end
        iter = iter + 1;     % increment minimization counter
        info.rp = rp * gam; % increase penalty multiplier (gamma = 5, here)
        U0 = Ustar;    % use current xstar as next x0
        info.Ts = Ustar(end);
%         if Ustar(end) < 0 
%             break;
%         end
        
    end

    %     [xi, yi, zi] = extractPos(xref, info.Nsteps+1);
    %     figure(k)
    %     plot3(xi,yi,-zi)
%     if Ustar(end) < 0
%         %must restart the loop with different parameters
%         info.Nsteps = info.Nsteps + 1;
%         info.Ts = 0.01;
%         Ustar = [Ustar(1:end-1); ones(info.dimM,1); info.Ts]; %create new guess for reshaped trajectory
%         redo_traj = true;
%     else
        x0 = info.x0;
        xref = predictStates(x0,Ustar,info);
        Gate_ID = [Gate_ID, k];
        Ts_total = [Ts_total, Ustar(end)];
        X_ref_total = [X_ref_total; xref];
        Ustar_total = [Ustar_total; Ustar(1:info.Nsteps*4)];
        Gatepos_total = [Gatepos_total, info.xf];
        toc/60
        redo_traj = false;
        k = k + 1;
%     end
    
    
    
    
    
end

Num_gates = gate.last;
Nstep = info.Nsteps;
Substep = info.substeps;

[xi, yi, zi] = extractPos(X_ref_total, length(X_ref_total)/12);
figure(k)
plot3(xi,yi,-zi)
hold on
scatter3(gate_x, gate_y, -gate_z)



























