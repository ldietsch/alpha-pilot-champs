function X = predictStates(x0, U, Ts, N, info)

X = zeros(12*(N+1),1);
m = 4;
k = 1;
j = 1;
l = 13;
t = 2;
X(1:12) = x0(1:12);
for i = 1:N
    
    xnew = forwardEuler(x0(l:t*12,1),U(j:k*m,1),Ts, info);
    X(l:12*t,1) = xnew;
    k = k + 1;
    j = j + m;
    l = l + 12;
    t = t + 1;
    x0(l:t*12) = xnew;
    
end

end