function [f] = costFunFixedTime(U, xref, info)

Nsteps = info.Nsteps;
m = info.dimM;


R = info.R;

Q = eye(length(xref));

f = 1/(Nsteps*info.Nstates)^2 *xref'*Q*xref + ...
    1/(info.dimM*Nsteps)^2 * U(1:Nsteps*m).'*R*U(1:Nsteps*m);