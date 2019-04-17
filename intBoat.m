function Y = intBoat(U, N)

x0 = [0;0];
X = x0;
Ts = U(1,N+1);

Y = [x0];
for i = 1:N-1
    u_step = U(:,i);
    dx = boat(u_step, Ts, i);
    X = X+dx;
    Y = [Y, X];
end