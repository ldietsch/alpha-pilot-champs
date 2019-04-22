function g = inequalConFixedTime(U, info)
Fupper = info.motorlimit;
Vupper = info.vel_upper;
gUpper = U./Fupper - 1;
gLower = -U;

g = [gUpper; gLower;];