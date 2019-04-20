function [g, h] = nlcons(U, x0, x_final_pos, quad_limits, info)

N = info.nMPC;
X = predictStates(x0, U, info);

h = (X((N)*12+1:N*12+3) - x_final_pos);

DX = predictDerivs(X, U, info);

v_max = quad_limits.v_max;
vel = extractVel(DX,N);
vel_cons = getVelCons(vel, v_max, N);
angles = extractAngles(X,N);
angleU = (angles - [pi;pi/2;pi])./[pi;pi/2;pi];
angleL = (-angles - [pi;pi/2;pi])./[pi;pi/2;pi];
g = [vel_cons angleU angleL];

end