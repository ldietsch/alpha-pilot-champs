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
input_file = '../input-files/test.xlsx';
[uMPC,~,exitflag,~] = ssnmpc(input_file,info);


%% Post-processing
info.post_processing = true;
x0 = zeros(n,1)*.2;
results = predictStates(x0, uMPC, info);
[x, y, z] = extractPos(results, info.Nsteps);
plot3(x,y,z,'--')

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


