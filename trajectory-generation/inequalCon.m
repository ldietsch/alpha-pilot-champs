function g = inequalCon(U, info)
Fupper = info.motorlimit;
Vupper = info.vel_upper;
Nsteps = info.Nsteps;
gUpper = U(1:end-1)./Fupper - 1;
gLower = -U(1:end-1);
Ts = U(end);
xREF = predictStates(info.x0,U,info);
% psi = extractHeading(xREF, Nsteps);
% psiUpper = (psi - 1)./pi;
% psiLower = (-psi - 1)./pi;

velF = xREF(end-9:end-7);
velFupper = velF./Vupper - 1;
velFlower = -velF./Vupper - 1;

% g = [gUpper; velFupper; psiUpper; Ts/0.05 - 1; gLower; velFlower; ...
%     psiLower; -Ts/0.0001 + 1];
g = [gUpper; velFupper; Ts/0.05 - 1; gLower; velFlower; ...
     -Ts/0.0001 + 1];