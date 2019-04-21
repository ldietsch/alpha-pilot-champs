function g = inequalConFixedTime(U, info)
Fupper = info.motorlimit;
Vupper = info.vel_upper;
Nsteps = info.Nsteps;
gUpper = U./Fupper - 1;
gLower = -U;
xREF = predictStates(info.x0,U,info);

velF = xREF(end-9:end-7);
velFupper = velF./Vupper - 1;
velFlower = -velF./Vupper - 1;

g = [gUpper; velFupper; gLower; velFlower];