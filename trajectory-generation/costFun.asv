function [f] = costFun(U, xref, info)

Nsteps = info.Nsteps;
m = info.dimM;
Ts = U(end);

R = info.R;

Q = eye(length(xref));

f = .01* xref'*Q*xref + .001 * U(1:Nsteps*m).'*R*U(1:Nsteps*m) + Nsteps*Ts^2;