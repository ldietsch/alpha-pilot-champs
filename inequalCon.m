function g = inequalCon(U, N)
Ts = U(1,N+1);

g = [];
for i = 1:N
g_i = [U(1)/7 - 1; -U(1);...
    U(2)/(pi) - 1; -U(2)/(pi) + 1];
g = [g, g_i];
end
g_ts = [Ts/0.1 - 1; -Ts/.001 + 1; 0;0];
g = [g, g_ts];
