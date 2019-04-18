function A = psuedoFun(U, info)
% info.Ts1 = U(end);

[h, xref] = equalCon(U, info);
f = costFun(U, xref, info);
g = inequalCon(U, info);


Nsteps = info.Nsteps;
lambda_g = info.lamg;
lambda_h = info.lamh;
rp = info.rp;

%fletchers principal

G = 0; H = 0;
for i = 1:Nsteps+2
    psi = max(g(i),-lambda_g(i)/2*rp);
    G = G + lambda_g(i)*psi + rp*psi^2;
end

for j = 1:3
    H = H + lambda_h(j)* h(j) + rp*h(j)^2;
end


A = f + G + H;