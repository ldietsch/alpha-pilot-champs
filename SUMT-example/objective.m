function f = objective(x0, U, info)
N = info.Nstates;
% Q = 
% R = 
X = predictStates(x0, U, info);
f = X'*Q*X + U(1:end-1)'*R*U(1:end-1) + N*U(end);


end
