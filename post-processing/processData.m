function X = processData(results_file,input_file)

x0 = zeros(12,1);
[~,~,data] = xlsread(results_file);
data = cell2mat(data(2:end,1:end));
U = data(1:end,2:end);
ID = data(1:end, 1);
X = zeros(size(U,1),12*size(U,2)/4);


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
[~,~,params] = xlsread(input_file,'params');

for i = 1:length(ID)
    %% reading in data from input file
    info.Nsteps = params{i+1,2};
    info.substeps = params{i+1,4};
    info.Ts = params{i+1,3};
    n = 12; %Number of states
    m = 4; %Number of inputs
    info.dimM = m;
    info.Nstates = n;
    info.rho = .1; %multiplier to control inputs in cost function
    info.post_processing = true;
    info.R = eye(info.dimM*info.Nsteps);
    X(i,:) = predictStates(x0,U(i,:)',info);
    
end

end

