function g = inequalCon(U, info)
Fupper = info.motorlimit;
Nsteps = info.Nsteps;
Ts = U(end);
% Ts = info.Ts1;
g = zeros(Nsteps*8+2,1);
j = 1;
for i = 1:1:Nsteps*4
    g(j,1) = U(i,1)/Fupper - 1;
    j = j+1;
    g(j,1) = -U(i,1);
    j = j+1;
end
g(end-1) = Ts/0.05 - 1;
g(end) = -Ts/0.001 + 1;