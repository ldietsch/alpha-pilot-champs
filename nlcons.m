function [g, h] = nlcons(U, X0, x_final_pos, Ts, N, quad_limits, info)

X = predictStates(X0, U, Ts, N, info);
h = (X((N-1)*12+1:(N-1)*12+4) - x_final_pos)./x_final_pos;

DX = predictDerivs(X, U, N, info);
a_max = quad_limits.a_max;
acc = extractAccel(DX,N);
accel_cons = getAccelCons(acc, a_max);

v_max = quad_limits.v_max;
vel = extractVel(DX,N);
vel_cons = getVelCons(vel, v_max);

g = [accel_cons; vel_cons];

end