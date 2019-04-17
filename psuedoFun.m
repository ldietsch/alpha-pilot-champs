function A = psuedoFun(U, rp, lambda_g, lambda_h, N)

f = costFun(U ,N);
g = inequalCon(U ,N);
X = intBoat(U, N);
h = equalCon(X, N);

%fletchers principal

G = 0; H = 0;
for j = 1:4
    for i = 1:N
        psi = max(g(j,i),-lambda_g(j,i)/2*rp);
        G = G + lambda_g(j,i)*psi + rp*psi^2;      
    end
end

for j = 1:2
H = H + lambda_h(j)* h(j) + rp*h(j)^2;
end


A = f + G + H;