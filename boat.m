function State = boat(u, Ts ,N)
w = 0.0008 * N^2;

State = [u(1)*cos(u(2))*Ts;...
    u(1)*sin(u(2))*Ts + w];