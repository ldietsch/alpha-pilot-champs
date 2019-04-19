function [uMPC, fval, exitflag, output, info] = ssnmpc(filename, info)
%ssnmpc is single-shooting nonlinear model predictive control
[~,~,input_file] = xlsread(filename);
num_gates = input_file{1,2};
segment_IDs = input_file{2,2:num_gates+1}; %will update all vectors according to actual file
Nsteps = input_file{3,2};
info.Nsteps = Nsteps;
info.substeps = input_file{4,2};
Ts = input_file{5,2:num_gates+1};
info.Ts = Ts(1);
Xref = cell2mat(input_file(6,2:3*num_gates*Nsteps+1))'; 
Uref = cell2mat(input_file(7,2:4*num_gates*Nsteps+1))'; 
gatePositions = cell2mat(input_file(8,2:3*num_gates+1))';
info.two_sample_times = false;

%nMPC = Nsteps for single-shooting
nMPC = Nsteps;
info.nMPC = nMPC; % the steps in the moving horizon
info.getting_next_state = false; %for one step only

Q0 = (1-info.rho).*eye(3);
Q = getQbar(Q0,3,nMPC);
R = eye(info.dimM*info.nMPC);
m = info.dimM;
A = []; b = []; Aeq = []; beq = []; 
a_max = 3*info.g;
lb = zeros(m*nMPC,1); 
m_tot = info.m;
ub = (m_tot*a_max/4)*ones(m*nMPC,1);
v_max = 80*1000/60^2; %80 km/h
quad_limits = struct('a_max',a_max,'v_max',v_max);

info.currentGate = 0;
n_total = num_gates; %total iterations is the number of segments for single-shooting
uMPC = zeros(numel(Uref),1);
options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',1e3);

j = 1;
k = 1;
l = 1;
    for i = 1:n_total
    info.iter = i;   
    info.nextGateID = floor(i/Nsteps) + 1;
    [cur_uref, cur_xref] = getNextTraj(Xref, Uref, info);    
    if (i + nMPC >= info.nextGateID*Nsteps)&&(info.nMPC~=info.Nsteps)
        info.two_sample_times = true;
        info.Ts1 = Ts(nextGateID);
        info.Ts2 = Ts(nextGateID+1);
%         x_pos = gatePositions(3*info.currentGate+1:3*info.nextGateID);   
        U0 = cur_uref;
        x0 = cur_xref(1:info.Nstates);
    else
        info.two_sample_times = false; %make sure it resets to false
%         x_pos = gatePositions(3*info.currentGate+1:3*info.nextGateID);
        x_pos = cur_xref(end-2:end);
        info.Ts1 = Ts(info.nextGateID);
        info.Ts2 = []; %reset to empty as a failsafe while debugging
        U0 = cur_uref;
        x0 = cur_xref(1:info.Nstates);        
    end
  
%     if (nMPC == n_total - i)
%         x_pos = gatePositions.n_segments;
%     elseif (nMPC > n_total - i) 
%         info.nMPC = n_total - i + 1;
%         nMPC = info.nMPC;
%         Q = getQbar(Q0,3,nMPC);
%         R = eye(info.dimM*info.nMPC);
%         lb = zeros(m*nMPC,1); 
%         ub = (m_tot*a_max/4)*ones(m*nMPC,1);
%     end
    
        %the minimization at each time step
        [U, fval, exitflag, output] = fmincon(@(U)costFunction(U, ...
            cur_uref, cur_xref, Q, R, x0, info),U0, A, b, Aeq, beq,...
            lb, ub, @(U)nlcons(U, x0, x_pos, quad_limits, ...
            info),options)
        
        uMPC(info.currentGate+1:m*info.Nsteps*info.nextGateID,1) = U';
        info.getting_next_state = true;        
        x0 = predictStates(x0, uMPC(l:l+3,1), info);
        info.getting_next_state = false;       
        info.currentGate = info.nextGateID;


        l = l + m;
        
    end   
    
end


