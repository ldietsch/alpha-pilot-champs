function [uMPC, fval, exitflag, output] = ssnmpc(filename, info)
%ssnmpc is single-shooting nonlinear model predictive control
[~,~,input_file] = xlsread(filename);
num_gates = input_file{30,2};
segment_IDs = input_file{2,1:num_gates}; %will update all vectors according to actual file
Nsteps = input_file{num_gates+2,2};
info.Nsteps = Nsteps;
info.substeps = input_file{num_gates+3,2};
Ts = input_file{:,3};
info.Ts = Ts(1);
Xref = input_file{2:num_gates+1,4:Nsteps*num_gates+4};
Uref = input_file{2:num_gates+1, Nsteps*num_gates+5:2*Nsteps*num_gates + 5};
gatePositions = input_file(num_gates+4:num_gates+6,1:num_gates); %should be a structure
info.two_sample_times = false;
nMPC = 10;
info.nMPC = nMPC; % the steps in the moving horizon

info.getting_next_state = false; %for one step only

Q0 = eye(3);
Q = getQbar(Q0,3,nMPC);
R = eye(info.dimM*info.nMPC);

A = []; b = []; Aeq = []; beq = []; 
a_max = 3*g;
lb = zeros(m*nMPC,1); 
ub = (m_tot*a_max/4)*ones(m*nMPC,1);
v_max = 80*1000/60^2; %80 km/h
quad_limits = struct('a_max',a_max,'v_max',v_max);

n_total = Nsteps*num_gates;
uMPC = zeros(numel(Uref));
options = optimoptions('fmincon','Algorithm','sqp');

j = 1;
k = 1;
l = 1;
    for i = 1:n_total
    info.iter = i;   
    info.nextGateID = floor(i/Nsteps) + 1;
    [cur_uref, cur_xref] = getNextTraj(Xref, Uref, info);    
    if i + nMPC >= nextGateID*Nsteps
        info.two_sample_times = true;
        info.Ts1 = Ts(nextGateID);
        info.Ts2 = Ts(nextGateID+1);
        x_pos = gatePositions.nextGateID;
    else
        info.two_sample_times = false; %make sure it resets to false
        info.Ts1 = Ts(nextGateID);
        info.Ts2 = []; %reset to empty as a failsafe while debugging
    end
  
    if (nMPC == n_total - i)
        x_pos = gatePositions.n_segments;
    elseif (nMPC > n_total - i) 
        info.nMPC = n_total - i + 1;
        nMPC = info.nMPC;
        Q = getQbar(Q0,n,nMPC);
        R = eye(info.dimM*info.nMPC);
        lb = zeros(m*nMPC,1); 
        ub = (m_tot*a_max/4)*ones(m*nMPC,1);
    end
    
        %the minimization at each time step
        [U, fval, exitflag, output] = fmincon(@(U)costFunction(U, ...
            cur_uref, cur_xref, Q, R, x0, info),U0, A, b, Aeq, beq,...
            lb, ub, @(U)nlcons(U, x0, x_pos, quad_limits, ...
            info),options);
        
        info.getting_next_state = true;
        uMPC(l:l+3) = U(1:m);
        x0 = predictStates(x0, uMPC(l:l+3), info);
        info.getting_next_state = false;
        l = l + m;
        
    end   
    
end


end

