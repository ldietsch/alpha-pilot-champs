function [f] = costFunFixedTime(U, xref, info)

Nsteps = info.Nsteps;
m = info.dimM;


R = info.R;
% Q0 = diag([1/4 1/4 1/4 1/2 1/2 1/2 1/4 1 1/4 1 1 1]);
Q = getQbar(Q0,info.Nstates,Nsteps);
% Q = eye(length(xref)-12);
% Q(end-11:end,end-11:end) = zeros(12);
xref = xref(13:end);
f = (1-info.rho).*xref'*Q*xref + ...
   info.rho * U'*R*U;