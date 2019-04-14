function X = predictStates(x0, U, info)

dimM = info.dimM;
substeps = info.substeps;

if info.post_processing %false by default
    N = info.Nsteps;
elseif info.optim_sol
    N = info.nMPC;
else
    N = 1;
end
X = zeros(12*(N+1),1);
xV = zeros(12, N*substeps+1);
xV(:,1) = x0;
X(1:12) = x0(1:12);
n = 1;
l = 13;

%this integration scheme follows the single shooting scheme
for i = 1:N
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