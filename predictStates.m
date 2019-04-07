function X = predictStates(x0, U, Ts, N, info)

X = zeros(12*N,1);
m = 4;
k = 1;
j = 1;
l = 1;
for i = 1:N
    
    xnew = forwardEuler(x0,U(j:k*m,1),Ts, info);
    X(l:12*i,1) = xnew;
    k = k + 1;
    j = j + m;
    l = l + 12;
    x0 = xnew;
    
end

end