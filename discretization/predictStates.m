function X = predictStates(x0, U, info)
%U is an nMPC*dimM x 1 vector and x0 is Nstates x 1 vector
dimM = info.dimM;
substeps = info.substeps;
if info.post_processing %false by default
    N = info.Nsteps;
elseif info.getting_next_state
    N = 1;
else
    N = info.nMPC;
end
X = zeros(12*(N+1),1);
xV = zeros(12, N*substeps);
xV(:,1) = x0;
X(1:12) = x0(1:12);
n = 1;
l = 13;

%this integration scheme follows the single shooting scheme
for i = 1:N
    if info.two_sample_times
        if i <= info.nextGateID*info.Nsteps - info.iter
            info.Ts = info.Ts1;
        else
            info.Ts = info.Ts2;
        end
    else
       info.Ts = info.Ts1; 
    end
        
    for j = 1:substeps
        k = (i - 1)*substeps + j;
        [xnew, ~] = RK45(xV(:,k),U(n:i*dimM,1), info);
        xV(:,k+1) = xnew;
    end
    n = n + dimM;
    X(l:l+11,1) = xnew;
    l = l + 12;
end

end