function [g, h] = nlcons(U, X0, x_final_pos, quad_limits, info)

N = info.nMPC;
X = predictStates(X0, U, info);
if ~isempty(x_final_pos)
    h = (X((N)*12+1:N*12+3) - x_final_pos);
else
    h =[];
end

DX = predictDerivs(X, U, info);

v_max = quad_limits.v_max;
vel = extractVel(DX,N);
vel_cons = getVelCons(vel, v_max, N);

g = vel_cons;

end