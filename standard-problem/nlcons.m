function [g, h] = nlcons(U, X0, x_final_pos, quad_limits, info)

N = info.nMPC;
X = predictStates(X0, U, info);

h = (X((N)*12+1:N*12+3) - x_final_pos);

DX = predictDerivs(X, U, info);

v_max = quad_limits.v_max;
vel = extractVel(DX,N);
vel_cons = getVelCons(vel, v_max, N);

g = vel_cons;

end