function [f] = MinFunc(X, U, Ts, N)
% Function to be minimized
rho = 1; Q_bar = eye(12);


f = X'*Q_bar*X + rho*U'*eye(4)*U + N*Ts;