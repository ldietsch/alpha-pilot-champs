function A = psuedoFunFixedTime(U, info)
% info.Ts1 = U(end);

[h, xref] = equalConFixedTime(U, info);
f = costFunFixedTime(U, xref, info);
g = inequalConFixedTime(U,xref, info);


Nsteps = info.Nsteps;
lambda_g = info.lamg;
lambda_h = info.lamh;
rp = info.rp;

%fletchers principal

G = 0; H = 0;
psi = max(g,-lambda_g/2*rp);
G = G + lambda_g'*psi + rp*(psi'*psi);



H = H + lambda_h'* h + rp*(h'*h);



A = f + G + H;