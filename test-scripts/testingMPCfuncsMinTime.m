clear all
close all
clc
%% adding all folders necessary for compilation
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
tol = 1e-4;
info = struct('m',m_tot,'g',g,'J',J,'l',l,'k_tau',k_tau,'tol',tol);
n = 12; %Number of states
m = 4; %Number of inputs
info.dimM = m;
info.Nstates = n;
info.nMPC = 10;
info.Nsteps = 40; % will be necessary for overall trajectory - total iterations will nMPC*Nsteps
info.substeps = 25; %steps to be used with integrator
info.rho = 1e-4; %multiplier to control inputs in cost function
nMPC = info.nMPC; %N steps ahead in MPC
info.optim_sol = true;
Ts = 0.01;
info.Ts = Ts; %time step to be used in integrator
info.post_processing = false;
Q0 = [diag([1 1 1])];
Q = getQbar(Q0,3,nMPC);
R = eye(info.dimM*info.nMPC);

x0 = zeros(n,1);
Uref = [ones(m*(nMPC),1)']';

Gates = readtable('Gate Locations.xlsx');
gate.N = Gates{:,1};
gate.x = Gates{:,2};
gate.y = Gates{:,3};
gate.z = Gates{:,4};
gate.current = 0;
gate.last = gate.N(end);

%% testing main MPC algorithm
info.post_processing = false;
A = []; b = []; Aeq = []; beq = []; 
a_max = 3*g;
lb = [zeros(m*nMPC,1); 1e-6]; 
ub = [(m_tot*a_max/4)*ones(m*nMPC,1); 0.05];
v_max = 80*1000/60^2; %80 km/h
quad_limits = struct('a_max',a_max,'v_max',v_max);


options = optimoptions('fmincon','Algorithm','interior-point',...
    'MaxFunctionEvaluations',1e4);

Gate_ID = [];
Ts_total = [];
X_ref_total = [];
Ustar_total = [];
Gatepos_total = [];
gate_x = []; gate_y = []; gate_z = [];
for i = 1:gate.last
    tic
    if i == 1
        info.x0 = zeros(12,1);        
    else
        info.x0 = xref(end-11:end);%previous states
    end
    nextGate = [gate.x(i+1),gate.y(i+1),gate.z(i+1)];
    info.xf = nextGate;
    
    U0 = ones(m*info.nMPC,1) * 2.6978;     %create initial input for minimization
    U0 = [U0; Ts];            %add Ts to input vector
    
    x0 = info.x0;
    
    info.optim_sol = true;
    [uMPC, fval, exitflag, output] = fmincon(@(U)costFunctionMinTime(U,Q,R,x0,...
    info),U0, A, b, Aeq, beq, lb, ub, ...
    @(U)nlcons(U, x0, info.xf, quad_limits, info),options);

    xref = predictStates(x0,uMPC,info);
    Gate_ID = [Gate_ID, i];
    Ts_total = [Ts_total, uMPC(end)];
    X_ref_total = [X_ref_total; xref];
    Ustar_total = [Ustar_total; uMPC(1:info.Nsteps*4)];
    Gatepos_total = [Gatepos_total, info.xf];
    toc/60
end

%% Post-processing
info.post_processing = true;
x0 = zeros(n,1)*.2;
% results = predictStates(x0, uMPC, info);
[x, y, z] = extractPos(X_ref_total, info.nMPC*gate.last);
plot3(x,y,-z,'--')
hold on
scatter3(gate_x, gate_y, -gate_z)
color_palette = {};
p_hand = plot3(x,y,z,'LineStyle','--','LineWidth',2);
hold on
color_palette = p_hand.Color;
plot3(x0(1),x0(2),x0(3),'^','MarkerSize',6,'MarkerEdgeColor',p_hand.Color,...
                                                 'MarkerFaceColor',p_hand.Color);
plot3(x_final_pos(1), x_final_pos(2),x_final_pos(3),'p','MarkerSize',8,'MarkerEdgeColor',p_hand.Color,...
                                                 'MarkerFaceColor',p_hand.Color);
           
plot3(x,y,z,'LineStyle','-','LineWidth',2,'Color',color_palette)
        hold on
xlim([(min(x)-0.5) (max(x)+0.5)])
ylim([(min(y)-0.5) (max(y)+0.5)])
zlim([(min(z)-0.5) (max(z)+0.5)])
start_hand = plot3(10000,10000,10000,'b^','MarkerSize',6);
end_hand   = plot3(10000,10000,10000,'bp','MarkerSize',8);
legend([start_hand,end_hand],{'Start Point','End Point'},'Location',...
    'northeast');
title("AlphaPilot Race Course")
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
    
drawnow
    
setFont();
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
grid ON
simTrajectories3D(x,y,z,color_palette);
