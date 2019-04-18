function [f] = costFun(U, xref, info)

Nsteps = info.Nsteps;
m = info.dimM;
Ts = U(end);

% Ts = info.Ts1;
R = info.R;

% x0 = info.x0;
% xref = predictStates(x0,U,info);

Q = eye(length(xref));

f = .001* xref'*Q*xref + .001 * U(1:Nsteps*m).'*R*U(1:Nsteps*m) + Nsteps*Ts^2;