function [g, h] = nlcons(U, X0, x_final_pos, Ts, N, quad_limits, info)

X = predictStates(X0, U, Ts, N, info);
if ~isempty(x_final_pos)
    h = (X((N)*12+1:N*12+3) - x_final_pos)./x_final_pos;
else
    h =[];
end

DX = predictDerivs(X, U, N, info);
a_max = quad_limits.a_max;
acc = extractAccel(DX,N);
accel_cons = getAccelCons(acc, a_max, N);

v_max = quad_limits.v_max;
vel = extractVel(DX,N);
vel_cons = getVelCons(vel, v_max, N);

g = [accel_cons; vel_cons];

end